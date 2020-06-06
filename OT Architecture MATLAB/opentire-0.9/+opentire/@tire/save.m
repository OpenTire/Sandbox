%% save - Save the tire object to a .TIR file

% part of mftire 1.1.0
function ok = save(tire, filename, header, format, compare_params)
  
  if nargin < 3, header = []; end
  if nargin < 4, format = 0; end
  if nargin < 5, compare_params = []; end
  
  ok = false;
  
  % Index values for special case parameter names
  PLUS = 2; % '+'
  TYRESIDE = 3;
  PROPERTY_FILE_FORMAT = 8;
  
  if (~tire.is_valid) %then
    fprintf('%s: Can''t save tire with invalid parameters\n', filename);
    return;
  end
  
  diff = false;
  
  info = tire.PARAM_INFO;
  
  defaults = info.defaults;
  if ~isempty(compare_params)
    diff = true;
  end
  
  if isempty(filename)
    fid = 1; % write to console
  else
    fid = fopen(filename, 'wt'); % text mode
    if fid < 0
      return;
    end
  end

  % print out header
  if ~diff
    if isempty(header)
      header = tire.mdi_header;
    end
    fprintf(fid, '%s', header);
  end
  
  diff_count = 0;
  tire_version = tire.model_version;
  
  for i=1:info.info_count
    version = info.versions(i);
    val_index = info.val_indexes(i);
    
    % [SECTION HEADER]
    if version == -1 && ~diff
      if format == 0 && i > 1
        fprintf(fid, '\n');
      end
      fprintf(fid, '%s\n', info.names{i});
      continue;
    end
    
    % values with negative indices don't exist
    if val_index < 1
      continue;
    end
    
    name = info.names{i};
    
    if val_index == PLUS
      value = tire.PLUS;
    elseif ~isprop(tire,name)
      % param_info table error: mftire_validate() needed
      fprintf('Parameter %s doesn''t exist\n', name);
      fprintf('param_info table error: run mftire_internal([], ''validate'');\n');
      return;
    else
      value = tire.(name);
    end
    
    if diff
      if value == compare_params(val_index)
        continue;
      end
      diff_count = diff_count + 1;
    else
      % don't write obsolete values if same as default
      if version < 0 && tire_version > abs(version)
        if value == defaults(val_index)
          continue;
        end
      end
    end

    % convert to printable value
    if val_index == TYRESIDE
      if (tire.TYRESIDE ~= 0)
        value = '''RIGHT''';
      else
        value = '''LEFT''';
      end
    elseif val_index == PROPERTY_FILE_FORMAT
      % Only 'USER' format supported by load/save
      value = '''USER''';
    else
      value = sprintf('%.8g', value);
    end
    
    if format == 0 % long format with comments
      name = [pad(name,7) ' = ' value];
      cmt = info.comments{i};
      
      if diff
        name = pad(name, 23);
        name = [name sprintf(' $ (%.8g)', compare_params(val_index))];
        name = pad(name, 39);
        if ~isempty(cmt)
          cmt = [' ' cmt];
        end
        diff_count = diff_count + 1;
      else
        name = pad(name, 23);
        cmt = sprintf(' $ %d %s', val_index, cmt);
      end
      fprintf(fid, '%s%s\n', name, cmt);
    elseif format == 1 % short format, no comments
      fprintf(fid, '%s = %s\n', name, value);
    end
  end
  
  if diff
    fprintf(fid, '\n$ tire.diff: %d parameters different\n', diff_count);
  end
  
  if fid > 2
    fclose(fid);
  end
  
  if diff
    ok = diff_count;
  else
    tire.filename = filename;
    ok = true;
  end
end
