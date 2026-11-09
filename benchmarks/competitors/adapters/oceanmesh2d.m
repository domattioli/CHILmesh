% Standalone Octave adapter for OceanMesh2D (GNU Octave attempt, no MATLAB seat assumed).
% Usage: octave --no-gui --eval "bench_mesh='in.14'; bench_out='out.txt'; run('oceanmesh2d.m')"
% Reads a fort.14 with OceanMesh2D's msh class, cleans it, writes points and elements as text.
% The OceanMesh2D call shape is unverified until research.md (T003/T005). Any error exits
% non-zero; the runner records `not-run` plus the stderr tail.
% OceanMesh2D must be on the Octave path (see README runbook).
try
  tic;
  m = msh(bench_mesh);            % OceanMesh2D mesh class
  m = clean(m, 'passive');        % OceanMesh2D mesh cleaning
  op_s = toc;
  dlmwrite([bench_out '.points'], m.p, 'precision', '%.17g');
  dlmwrite([bench_out '.cells'], m.t, 'precision', '%d');
  printf('{"tool":"OceanMesh2D","op_s":%.6f}\n', op_s);
catch err
  fprintf(stderr, '%s\n', err.message);
  exit(1);
end_try_catch
