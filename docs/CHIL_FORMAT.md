# The `.chil` container format

**EXPERIMENTAL.** The `.chil` format and the `chilmesh.chil_io` module are experimental. The current layout is `format_version = 0`. CHILmesh 1.8.x reads `format_version = 0` archives written by 1.8.x; later minor releases may require migration.

## What a `.chil` is

A `.chil` file is a zip archive. Its first entry is `manifest.json`. Every other entry is a part listed in that manifest. A `fort.14` is always present. The manifest records a size and a SHA-256 for every part, so a reader can tell whether any byte changed after packing. There are no directory entries, no encryption, and the zip comment and trailing bytes are ignored on read and empty on write. Written entries use DEFLATE; STORED and DEFLATE are accepted on read.

The `.chil` suffix is matched case-insensitively, so `.CHIL` works for save, load and summary. The `validate` command and `validate_chil` do not look at the suffix.

## Manifest

`manifest.json` is UTF-8 JSON, one object, no BOM. The writer emits sorted keys, 2-space indent, ASCII escapes, a trailing newline, and refuses NaN and infinity. Any top-level key not in this table is refused (`unknown-key`).

| Key | Type | Required | Rule |
|---|---|---|---|
| `format_version` | integer | yes | exactly `0` |
| `crs` | string | yes | 1 to 256 characters, no control characters; default `"unknown"` |
| `parts` | array | yes | 1 to `max_members - 1` entries, names unique |
| `domain_boundary` | string | no | the `name` of the part whose kind is `domain_boundary` |
| `summary` | object | no | cached `n_nodes`, `n_elements`, `n_triangles`, `n_quads`, `grid_name`, and optional `bbox` |
| `source_archive_sha256` | string | no | 64 lowercase hex digits; written only by `rewrite_chil` |
| `ext` | object of objects | no | extension tables, one per namespace |

Each part has exactly the keys `name`, `role`, `kind`, `size` and `sha256`. `size` is the number of decompressed bytes. `sha256` is the digest of the decompressed bytes.

The `crs` value is a label. CHILmesh stores it and returns it. It never interprets it and never transforms coordinates.

## Part kinds, roles and member names

| kind | member name | allowed roles |
|---|---|---|
| `fort14` | `fort.14` | `source`, `generated` |
| `fort13` | `fort.13` | `source` |
| `fort15` | `fort.15` | `source` |
| `domain_boundary` | `domain_boundary/<basename of the supplied file>` | `source` |
| `ext_member` | `ext/<owner>/<rest>` | `extension` |

Exactly one `fort14` part is required. At most one part each of `fort13`, `fort15` and `domain_boundary` is allowed. Parts of kind `ext_member` are limited only by `max_members`.

Roles describe how the bytes were made:

- `source`: bytes supplied by the caller and copied unchanged (made by `pack`).
- `generated`: bytes the `fort.14` writer produced from an in-memory mesh (made by `save`).
- `extension`: opaque bytes supplied by an extension owner. CHILmesh lists and hashes them and never reads them. Only `ext_member` parts have this role.

A role is a claim about how the bytes were made, not about who made them. If any part is `generated`, every other part must be an `ext_member`.

An extension owner adds a table `ext.<owner>` to the manifest and its own members under `ext/<owner>/`. An `ext_member` whose owner has no table is refused. Owner names match `^[a-z][a-z0-9_]{0,63}$` and may not equal a core manifest key. `save` and `pack` never write extension tables or members; only `rewrite_chil` does.

Member names are refused (`member-name`) when they are empty, absolute, contain a backslash, NUL or control character, start with a drive letter, have an empty, `.` or `..` segment, or exceed 255 UTF-8 bytes. Two names that are equal after Unicode case folding collide (`duplicate-name`).

## Pack and save

Both create an archive. They differ in how the `fort.14` bytes come to exist.

- `pack` is a photocopy. It reads existing files in binary mode and stores their bytes unchanged, CRLF included. All parts get role `source`. `pack` accepts `fort13`, `fort15` and a Domain Boundary file as companions.
- `save` is a retyped page. It writes the in-memory mesh through the `fort.14` writer, so the stored text is what CHILmesh writes, not what the original file contained. The `fort.14` part gets role `generated`. A loaded mesh keeps no reference to the archive it came from, so a saved archive never matches the loaded file byte for byte, and companions, extension tables and `source_archive_sha256` are dropped.

Byte-identical saves hold only within one Python and zlib build. Two saves of the same mesh on different builds may differ in compressed bytes.

## Load

`CHILmesh.load("x.chil")` validates the archive, decodes the `fort.14` as UTF-8, and parses it. A parse failure becomes `ChilFormatError` with check `fort14-parse`. A `.chil` carries triangle and quad meshes only: an element row with a node count other than 3 or 4 is refused as `fort14-parse` (naming the element id) in `pack`, validation, load and the summaries. Load writes nothing to disk. Companions and extension members are not parsed; read them with `read_member`.

## Summary

`summary_chil` and `chilmesh.summary` return a dict. The CLI `chilmesh summary` prints it.

- Plain mode opens only `manifest.json`. Cached values (`n_nodes`, `bbox` and so on) are declared by the manifest and not checked against the payload. The result has `verified = False`.
- Deep mode (`deep=True`, CLI `--deep`) validates the whole archive, parses the `fort.14`, and recomputes the counts and bounding box. The result has `verified = True`, a `recomputed_summary`, and a sorted `mismatches` list naming each cached field that differs. A `fort.14` with a non-finite coordinate is refused in deep mode as `fort14-parse`.

## Rewrite

`rewrite_chil(src, dst, ...)` is the only write path that keeps extension tables, extension members and untouched members. It validates `src`, copies every part not named in the request while re-checking its digest, applies the requested changes, sets `source_archive_sha256` to the SHA-256 of the source archive file, and validates the result before writing. `dst` may equal `src`.

It can set `ext` tables, add or replace `fort13`, `fort15`, `domain_boundary` and ext members, remove those parts, and change `crs`. It never replaces or removes the `fort.14`; such a request raises `request` and writes nothing. The reason: extension tables, signatures and companions describe one specific mesh, and swapping the mesh would leave them describing a mesh that is no longer in the archive.

## Validate

`chilmesh validate FILE` and `validate_chil(path)` check the whole archive and either return a report or raise `ChilFormatError`. The `.check` attribute holds one of these slugs.

| Slug | Raised for |
|---|---|
| `not-zip` | not a zip, bad central directory, CRC failure, local header mismatch |
| `limit-members`, `limit-member-bytes`, `limit-total-bytes`, `limit-manifest-bytes` | a size or count limit exceeded |
| `compression`, `encrypted` | unsupported entry compression or encryption |
| `member-name`, `duplicate-name` | unsafe or colliding member name |
| `manifest-missing`, `manifest-not-first` | no `manifest.json`, or not the first entry |
| `manifest-json`, `duplicate-key`, `non-finite` | invalid JSON, repeated key, NaN or infinity |
| `prototype-container` | top-level key `c1_version` from an earlier prototype |
| `format-version` | `format_version` is not `0` |
| `unknown-key` | a manifest, part or summary key not in the format |
| `manifest-schema` | wrong type or value for a known field |
| `part-role`, `part-kind`, `part-count` | invalid role or kind, bad role and kind pairing, wrong number of parts |
| `ext-member-name`, `ext-member-owner` | bad `ext/<owner>/<rest>` name, or owner without a table |
| `unlisted-member`, `missing-member` | a zip entry not in the manifest, or a listed part not in the zip |
| `domain-boundary` | `domain_boundary` key and part do not agree |
| `size-mismatch`, `digest-mismatch` | decompressed bytes differ from the manifest size or SHA-256 |
| `node-id-zero` | node id 0 in an element row or a boundary row |
| `node-id-range` | node id below 1 or above `n_nodes` |
| `archive-changed` | the archive file changed while it was being validated |
| `fort14-parse` | the `fort.14` cannot be parsed (pack, load, deep summary), or an element row has a node count other than 3 or 4; the message names the element id |
| `ext-namespace`, `ext-table` | bad extension namespace, or a table that is not an object |
| `request` | an invalid write request (bad CRS label, a `bytes` path, bad kind, changing the `fort.14` in a rewrite, input not a regular file, input changed between passes) |

Pack and validate accept any `fort.14` whose header, elements and node ids pass their checks, so a packed archive can validate and still be refused by load or `--deep` for a reason those checks do not cover.

Node id 0 is refused everywhere. A `fort.14` is one-based; id 0 in an element or boundary row raises `node-id-zero`. The boundary scan reads rows in the reader's order and checks every id it reads up to the first row that does not parse; rows after that row are not scanned.

Memory use is bounded: the boundary rows of a `fort.14` are scanned as they stream in and are not buffered, so a boundary block of any size is scanned in full.

Extension tables are preserved by value in the JSON data model: the same keys, types and values, with floats as IEEE doubles (for example `1E2` reads back equal to `100.0`). They are not preserved as text. Content that must stay byte for byte belongs in an extension member (`ext/<owner>/...`), which is kept byte for byte.

`summary.grid_name` is the first line of the `fort.14` (ended by `\n`, `\r\n` or `\r`), stripped, cut to 256 UTF-8 bytes at a character boundary.

The manifest must be the first entry in the archive bytes (local header at offset 0), not only first in the central directory; otherwise `manifest-not-first`. A Domain Boundary basename is at most 200 UTF-8 bytes (`member-name`). Paths given as `bytes` are refused with `request`. A `crs` given to a writer must be 1 to 256 characters without control characters, else `request`.

## Limits and defaults

`ChilLimits` is a frozen dataclass. Every function that reads or writes an archive takes `limits=`. Limits apply on read and on write; a write that would produce an archive over a limit raises the same `limit-*` check before any output file exists.

| Field | Default | Counts |
|---|---|---|
| `max_members` | 64 | zip entries, manifest included |
| `max_member_bytes` | 2 GiB | decompressed bytes per member, manifest excluded |
| `max_total_bytes` | 4 GiB | decompressed bytes over all members, manifest included |
| `max_manifest_bytes` | 16 MiB | decompressed manifest bytes |

## Concurrency

Writes are not safe to run from several threads of one process at the same moment.

## Out of scope for this version

These belong to later phases of the `.chil` roadmap (design thread: issue #201):

- Valence fingerprints, registry lookup and the meaning of signatures: phase P2. CHILmesh carries Valence's `ext.valence` table and `ext/valence/` members as opaque bytes and does not interpret them.
- Project tie-in such as "produced by" sections and boundary outlines: phase P3.
- UGRID export and a standalone reader: phase P4.
- Binary array storage for coordinates and connectivity: phase P5.
- CRS transforms: not part of format version 0. The `crs` field is a label only.
