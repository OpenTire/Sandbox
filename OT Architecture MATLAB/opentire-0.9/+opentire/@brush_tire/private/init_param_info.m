%% init_param_info - Initialize the param_info structure
% 
% part of mftire 1.1.0
function [info, defs] = init_param_info()
% First get param_defs, a cell array in the following format:
%
% Column 1: Parameter index (1-based)
% Column 2: MF model version 520, 600, 610, etc.
%           Negative means the parameter was obsolete in that version
%           (and will not be included or saved in later model versions)
% Column 3: Default value.  Usually 0, sometimes 1, but can be any value.
% Column 4: Parameter name.
% Column 5: Parameter description/comment.

  defs = param_defs();
  prev_param_count = 279;
  prev_table_size = 334;
  
  % Renumber parameter indices in first column,
  % and count number of parameters >= 0.
  param_count = 0;
  table_size = size(defs,1);
  changed = (table_size ~= prev_table_size);
  for i=1:table_size
    if defs{i,1} >= 0
      if defs{i,1} ~= param_count + 1
        changed = true;
      end
      defs{i,1} = param_count + 1; % 1-based indices
      param_count = param_count + 1;
    end
  end

  if nargout == 2
    info = [];
    return;
  end
  
  % Print out new table if something changed.
  % Copy the text in the Command Window and paste it into
  % this file above.
  if param_count ~= prev_param_count || changed
    dump_param_info(defs);
    fprintf('\nparam_info parameter table has changed: ');
    fprintf('param_count = %d, was %d\n', param_count, prev_param_count);
    fprintf('table_size = %d, was %d\n', table_size, prev_table_size);
    fprintf('then change prev_param_count and prev_table_size lines in init_param_info.m\n');
  end
  
  info.info_count = size(defs,1);
  info.param_count = param_count;
  info.val_indexes = [defs{:,1}]';
  info.versions = [defs{:,2}]';

  % defaults table is indexed by val_index, not param_info row index
  vi = info.val_indexes;
  info.defaults = [defs{:,3}]';
  info.defaults = info.defaults(vi>0);
  
  info.names = upper(defs(:,4));
  info.comments = defs(:,5);
  
  % map from name to info table index
  info.name_map = containers.Map(info.names, 1:info.info_count);
  
  % map from value index to info index
  val_index_map = zeros(param_count,1);
  j=1;
  for i=1:size(defs,1)
    if defs{i,1} >= 0
      val_index_map(j) = i;
      j = j + 1;
    end
  end
  info.val_index_map = val_index_map;
  
end % init_param_info

function param_count = dump_param_info(defs)
  
  fprintf('\n\n  defs = {\n');
  param_count = 0;
  for i=1:size(defs,1)
    val_index = defs{i,1};
    model = defs{i,2};
    if val_index > 0 && model > 0
      param_count = param_count + 1;
    end
    
    default_value = defs{i,3};
    name = ['''' defs{i,4} ''''];
    comment = defs{i,5};
    fprintf('    %3d %4d %g %s ''%s''\n', val_index, model, default_value, pad(name, 8), comment);
  end
  fprintf('  };\n\n');
  fprintf('Copy new table text from the Command Window and paste into @tire/private/param_defs.m file\n');
end

