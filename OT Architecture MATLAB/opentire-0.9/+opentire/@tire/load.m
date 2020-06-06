%% load - Load a tire object from a .TIR file

% part of mftire 1.1.0
function [ok, header, params] = load(tire, filename)
  
  tire.filename = filename;
  tire.is_valid = false;
  
  ok = false;
  header = [];
  
  % Index values for special case parameter names
  PLUS = 2; % '+'
  TYRESIDE = 3;
  PROPERTY_FILE_FORMAT = 8;
  
  info = tire.PARAM_INFO;
  
  % Initialize with default values from param_info table
  params = info.defaults;
  for i=1:length(info.names)
    val_index = info.val_indexes(i);
    if val_index == PLUS
      tire.PLUS = params(val_index);
    elseif val_index == PROPERTY_FILE_FORMAT
      tire.PROPERTY_FILE_FORMAT = params(val_index);
    elseif val_index > 0
      tire.(info.names{i}) = params(val_index);
    end
  end
  
  if isempty(filename)
    return;
  end
  
  fid = fopen(filename, 'rt'); % text mode
  if fid < 0
    %warning('mftire: Can''t open tire file %s', filename);
    fprintf('mftire: Can''t open tire file %s\n', filename);
    return;
  end
  
  header = '';     % empty header
  current_section = '';   % waiting for first section
  skip_section = false;
  
  first_section = info.names{1};   % should be '[MODEL]'
  while ~feof(fid)
    line = fgetl(fid);
    if line == -1
      break;
    end
    line = strtrim(line);
    
    if isempty(line) || line(1) == '!' || line(1) == '$'
      if isempty(current_section)
        header = [header line newline];
      end
      continue;
    end
    
    if isempty(current_section)
      if contains(upper(line), first_section)
        current_section = first_section;
      else
        header = [header line newline];
      end
      continue;
    end
    
    if line(1) == '['
      end_index = strfind(line, ']');
      if isempty(end_index)
        end_index = length(line);
      end
      name = upper(line(1:end_index(1)));
      current_section = name;
      
      skip_section = false;
      if ~info.name_map.isKey(name)
        fprintf('%s: Section %s not supported\n', filename, name);
        % clear current section to skip everything in it
        skip_section = true;
      end
      
      %fprintf('section %s\n', name);
      continue;
    end
    
    if (skip_section)
      continue;
    end
    
    name_index = strfind(line, '=');
    if isempty(name_index)
      fprintf('%s: Parameter error: missing ''='': ''%s''\n', filename, line);
      continue;
    end
    name = upper(strtrim(line(1:name_index(1)-1)));
    
    value_index = strfind(line, '$');
    if isempty(value_index)
      value_index = length(line)+1;
    end
    
    value = strtrim(line(name_index(1)+1:value_index(1)-1));
    
    % strip quotes from string values
    if value(1) == '''' && value(end) == ''''
      value = value(2:end-1);
    end
    
    %fprintf('%s = %s\n', name, value);
    if ~info.name_map.isKey(name)
      fprintf('%s: Parameter %s = %s doesn''t exist\n', filename, name, value);
      continue;
    end
    
    index = info.val_indexes(info.name_map(name));
    if index < 1
      fprintf('%s: Parameter %s = %s not supported\n', filename, name, value);
      continue;
    end
    
    if index == PLUS
      % Assign to special case property PLUS instead of '+'
      value = str2double(value);
      tire.PLUS = value;
    elseif index == TYRESIDE
      % Convert TYRESIDE string value to a boolean (true == RIGHT)
      if (isempty(value) || upper(value(1)) == 'L')
        value = 0;
      else
        value = 1;
      end
      tire.TYRESIDE = value;
    elseif index == PROPERTY_FILE_FORMAT
      % Check for the only PROPERTY_FILE_FORMAT we support
      if (~strcmpi(value,'USER'))
        fprintf('%s: PROPERTY_FILE_FORMAT %s not supported, using ''USER''\n', filename, value);
      end
      value = 0;
      tire.PROPERTY_FILE_FORMAT = value;
    else
      value = str2double(value);
      if ~isfinite(value)
        fprintf('%s: Parameter value syntax error: ''%s''\n', filename, line);
      end
      if ~isprop(tire,name)
        fprintf('%s: Property %s doesn''t exist: fix properties in tire.m\n', filename, name);
      else
        tire.(name) = value;
      end
    end
    
    params(index) = value;
  end
  
  fclose(fid);
  
  % recompute internal variables
  
  tire.is_valid = true;
  tire.on_params_changed();
  
  ok = true;
end

