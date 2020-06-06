% utility - misc internal mftire functions for software development, etc

% part of mftire 1.1.0
function result = utility(tire, command)
  
  if nargin == 0
    validate();
    return;
  end
  
  switch lower(command)
    case 'dump_matlab_props'
      % write out matlab properties for pasting into opentire.mftire.m
      dump_matlab_props(tire, true)
      dump_matlab_props(tire, false);
      result = true;
      
    case 'validate'
      result = validate();
      
    otherwise
      fprintf('Unsupported utility() command: ''%s''\n', command);
      result = false;
  end
end

% write out all parameters and their values in
% matlab code format, for pasting into tire.m.
function param_count = dump_matlab_props(tire, cpp_format)
  
  info = tire.PARAM_INFO;
  
  defaults = info.defaults;
  
  for i=1:info.info_count
    version = info.versions(i);
    val_index = info.val_indexes(i);
    
    if version == -1 % section header
      if i > 1
        fprintf('\n');
      end
      if cpp_format
        fprintf('  // %s\n', info.names{i});
      else
        fprintf('%%%s\n', info.names{i});
      end
      continue;
    end
    
    % don't save values with negative indices
    if val_index < 1
      continue;
    end
    
    name = info.names{i};
    
    if ~isprop(tire,name)
      %fprintf('NOT DEFINED %s\n', name);
      %assert(false);
      %value = values(val_index);
      value = defaults(val_index);
    else
      value = tire.(name);
    end
    
    value = sprintf('%.7g', value);
    
    name = [pad(name,7) ' = ' value];
    
    if cpp_format
      name = [name ';'];
      cmt = info.comments{i};
      if ~isempty(cmt)
        if length(name) < 21
          name = pad(name, 21);
        end
        cmt = sprintf(' // %d %s',  val_index, cmt);
      end
      fprintf('  T %s%s\n', name, cmt);
    else
      cmt = info.comments{i};
      if ~isempty(cmt)
        if length(name) < 23
          name = pad(name, 23);
        end
        cmt = sprintf(' ;%% %d %s',  val_index, cmt);
      end
      fprintf('%s%s\n', name, cmt);
    end
  end
  
  param_count = info.param_count;
end

function ok = validate()
  
  tire = opentire.mftire(0);
  %print_dependent_props(tire);
  %print_setter_getters(tire);
  %print_rep_strings(tire);
  %return;
  
  % make sure loading, saving, then loading again doesn't lose values
  % (must be run from mfeval parent)
  file1 = 'tires/reference_tire_tvd3.tir';
  file1 = 'tires/MagicFormula62_Parameters.tir';
  
  [ok, header] = tire.load(file1);
  params = tire.get_param_vector();
  file2 = 'temp.tir';
  
  % write something into the header
  header = [header '! file: ''' file2 '''' char([13 10])];
  
  tire.save(file2, header);
  [ok2, header2] = tire.load(file2);
  params2 = tire.get_param_vector();
  
  if norm(params2 - params) ~= 0
    fprintf('Parameters differ after save/load round trip\nindices: \n');
    bad = find(params2~=params);
    fprintf('%d ', bad);
    fprintf('\n');
    ok = false;
  else
    fprintf('param_info table validation successful\n');
    ok = true;
  end
  
  delete('temp.tir');
  
end

% uncomment line above and call mftire_internal() with no parameters to print switch statements
function print_setter_getters(tire)
  info = tire.PARAM_INFO;
  val_name_map = info.val_index_map;
  names = info.names;
  fprintf('  for i=1:length(indices)\n');
  fprintf('    value = values(i);\n');
  fprintf('    switch indices(i)\n');
  for i=1:info.param_count
    fprintf('    case %d\n', i);
    fprintf('        t.%s = value;\n', names{val_name_map(i)});
  end
  fprintf('    otherwise\n');
  fprintf('    end\n');
  fprintf('  end\n');
  
  fprintf('\n\n');
  
  fprintf('  for i=1:length(indices)\n');
  fprintf('    switch indices(i)\n');
  for i=1:info.param_count
    fprintf('    case %d\n', i);
    fprintf('        value = t.%s;\n', names{val_name_map(i)});
  end
  fprintf('    otherwise\n');
  fprintf('    end\n');
  fprintf('    values(i) = value;\n');
  fprintf('  end\n');
end

function print_dependent_props(tire)
  info = tire.PARAM_INFO;
  val_name_map = info.val_index_map;
  names = info.names;
  
  for index=1:info.param_count
    name = names{val_name_map(index)};
    fprintf('function set.%s(t,p), t.P(%d) = p; end\n', name, index);
    fprintf('function p = get.%s(t), p = t.P(%d); end\n', name, index);
  end
end

% dump strings that can be used to replace t.NAME with tp(123) in .m files
function print_rep_strings(tire)
  info = tire.PARAM_INFO;
  names = info.names;
  val_indices = info.val_indexes;
  rep = cell(info.param_count, 2);
  i = 1;
  for ii=1:length(names)
    if val_indices(ii) <= 0
      continue;
    end
    rep{i,1} = names{ii};
    rep{i,2} = val_indices(ii);
    i = i + 1;
  end
  [~,ii] = sort(rep(:,1));
  ii = ii(end:-1:1);
  rep = rep(ii,:);
  for i=1:length(ii)
    fprintf('rep(''t.%s'',''tp(%d)'');\n',rep{i,1},rep{i,2});
  end
end
