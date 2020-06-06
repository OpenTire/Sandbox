%% param_defs - MF .TIR file parameter definition table.
%
% Edit this file to add or modify parameter definitions.
% You can add new parameters anywhere in the table,
% BUT the parameter index must always be the highest
% index in the table plus 1 -- this ensures that
% parameter indexes for existing parameters never change.
%
% Column 1: Parameter index (1-based)
% Column 2: MF model version 520, 600, 610, etc.
%           Negative means the parameter was obsolete in that version
%           (and will not be included or saved in later model versions unless non-zero)
% Column 3: Default value.  Usually 0, sometimes 1, but can be any value.
% Column 4: Parameter name.
% Column 5: Parameter description/comment.
%
% part of mftire 1.1.0
function defs = param_defs()
  defs = {
     -1   -1 0 '[MODEL]' ''
      1  520 0 'FITTYP' 'Magic Formula version number'
      2  610 0 '+'      'Moment representation 0=ground frame 1=wheel frame (not supported)'
      3  520 0 'TYRESIDE' 'Position of tire during measurements'
      4  520 0 'LONGVL' 'Reference speed'
  };
end
