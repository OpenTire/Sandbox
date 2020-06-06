% mftire_mfeval - Compare mftire output with MFeval 4.0 output.
% If called with no arguments, compare multiple tires from 'tires' folder.

% part of mftire 1.1.0
function mftire_mfeval(tire)
  if nargin == 0
    % run comparison on files in 'tires' folder
    compare_all();
    return;
  end
  
  want_turn_slip = false;
  
  rng 'default'; % reset random number generator for repeatable runs
  
  timing = true; % enable timing test
  
  compare_tol = 1e-4; % relative error compare tolerance

  % number of samples to compare
  samples = 10^3;
  %samples = 5;
  
  % list of timing sample counts to use
  % (set samples to 10^5 above and comment out 
  % '[1 samples]' line below to see how
  % AMAZINGLY FAST MFeval is when evaluating
  % 100,000 samples at once!
  timing_counts = 10.^[0 .5 1 1.5 2 2.5 3 3.5 4 4.5 5];
  % samples = 10^5
  
  % Just the extremes, no plot.
  timing_counts = [1 samples];
  
  if ischar(tire)
    tire = opentire.mftire(tire);
  end
  if ~tire.is_valid()
    fprintf('Invalid tire parameters\n');
    return; % file not found or other error
  end
  tire_file = tire.filename();
    
  close all;
  
  fprintf('\ntire file: %s\n', tire_file);
  
  
  % Load MFeval tire parameters using modified version
  % of MFeval 1.5 loadTIR.m
  mfe_t = mftire_loadTIR(tire_file);
  
  % If first turnslip coefficient is non-zero, use it.
  if (want_turn_slip)
    want_turn_slip = (tire.PDXP1 ~= 0);
  end
  fprintf('want_turn_slip = %d\n', want_turn_slip);
  
  % patch bad fields
  tire = validate_parameters(tire, tire_file);
  mfe_t = validate_parameters(mfe_t, tire_file);
  
  % Generate n random input samples
  Fz       = randrange(tire.FNOMIN/4, tire.FNOMIN*3, samples);
  kappa    = randrange(-1, 1, samples);
  alpha    = randrange(-20*pi/180,20*pi/180,samples);
  gamma    = randrange(-6*pi/180,6*pi/180,samples);
  phi_t    = randrange(1/-10,1/10,samples);
  kmh_to_ms = 1000/60/60;
  Vcx      = randrange(5*kmh_to_ms,150*kmh_to_ms,samples);
  psi_to_pa = 6895;
  
  press_range = [tire.PRESMIN tire.PRESMAX];
  % If no pressure info in file, use constant pressure of 0.
  if (press_range(1) >= press_range(2) || tire.NOMPRES == 0)
    press = zeros(1,samples);
  else
    % if pressure range is crazy, limit it to +/- some range
    press_delta = 0.20;
    press_range(1) = max(press_range(1), tire.NOMPRES * (1-press_delta));
    press_range(2) = min(press_range(2), tire.NOMPRES * (1+press_delta));
    
    press = randrange(press_range(1), press_range(2), samples);
  end
  omega = ones(1,samples)*nan;   % to be computed
  
  options = tire.new_options();
  options.want_turn_slip = want_turn_slip;
  
  USE_MODE = 211; % no limit checks, combined, no turn slip
  if options.want_turn_slip
    USE_MODE = USE_MODE + 1;    % include turn slip
  end
  
  INPUTS = [Fz(:) kappa(:) alpha(:) gamma(:) phi_t(:) Vcx(:) press(:) omega(:)];
  
  INPUTS = INPUTS(:,1:end-1); % remove omega
  
  OUTPUT = tire.mfeval(INPUTS, options);
  OUTPUT2 = mfeval(mfe_t, INPUTS, USE_MODE);
  
  % We only compare the first 20 MFeval 1.5-compatible outputs
  ncol = 20; % mftire 1.5 compatible
  
  OUTPUT = OUTPUT(:,1:ncol); % truncate extra outputs
  OUTPUT2 = OUTPUT2(:,1:ncol);
  
  col_scale = max( max(abs(OUTPUT)), max(abs(OUTPUT2)) );
  col_scale = max(col_scale, 1);
  
  d = abs(OUTPUT - OUTPUT2);
  
  % compute a relative error
  d = d ./ col_scale;
  %d = d./max(abs(OUTPUT),eps);
  
  % find the input row resulting in maximum error
  worst = max(max(d));
  i_worst = find(d(:)==worst);
  i_worst = i_worst(1); % in case of dupes
  worst_col = floor((i_worst-1)/samples) + 1;
  worst_row = (i_worst-1) - (worst_col-1)*samples + 1;
  
  worst_input = INPUTS(worst_row,:);
  
  avg = mean(d(:));
  d = max(d);
  
  % print any significant differences
  fprintf('\n%d samples compared: max abs( (mftire-MFeval)/MFeval ) = %.3g, average %.3g, worst row %d, col %d\n', samples, worst, avg, worst_row, worst_col);
  if (~isfinite(worst) || worst > compare_tol)
    
    fprintf('worst output (row %d, col %d):\n', worst_row, worst_col);
    fprintf('%10.6g ', OUTPUT(worst_row,:));
    fprintf('\n');
    fprintf('%10.6g ', OUTPUT2(worst_row,:));
    fprintf('\n\n');
    
    fprintf('worst input:\n');
    fprintf('%10.6g ', worst_input);
    fprintf('\n\n');
    
    e = sqrt(sum(d.^2,1)./size(d,1));
    
    fprintf('worst value (row %d, col %d): %.8g %.8g (rel err %.2g)\n\n', ...
      worst_row, worst_col, OUTPUT(worst_row, worst_col), OUTPUT2(worst_row, worst_col), e(worst_col));
    
    fprintf('rel error: ');
    fprintf('%8.2g', e);
    fprintf('\n');
    
    % loop so we can set breakpoints in both functions
    % to debug differences in results
    
    for i=1:10
      out = tire.mfeval(worst_input, options);
      out_mfe = mfeval(mfe_t, worst_input, USE_MODE);
    end
  else
    fprintf('**** COMPARISON OK\n');
  end
  
  if ~timing
    return;
  end

  % Timing test
  
  timing_counts = timing_counts(timing_counts <= samples); % remove counts we don't have samples for
  
  times = zeros(2,size(timing_counts,1));
  for ic = 1:length(timing_counts)
    timing_samples = round(min(samples, timing_counts(ic)));
    niter = round(samples / timing_samples);
    niter = max(niter,10);
    
    TINPUTS = INPUTS(1:timing_samples,:);
    
    fprintf('\nTiming: %d iterations of %d samples\n', niter, timing_samples);
    
    tic;
    for iter=1:niter
      OUTPUT = tire.mfeval(TINPUTS, options);
    end
    dt = toc*1e6/niter/timing_samples;
    times(ic,1) = dt;
    
    fprintf('mftire: %6.2f usec per eval\n', dt);
    
    tic;
    for iter=1:niter
      OUTPUT = mfeval(mfe_t, TINPUTS, USE_MODE);
    end
    dt = toc*1e6/niter/timing_samples;
    times(ic,2) = dt;
    fprintf('MFeval: %6.2f usec per eval\n', dt);
  end
  
  if length(timing_counts) > 2
    % plot eval times vs. log10(counts)
    close all;
    figure; hold on; grid on;
    plot(log10(timing_counts), log10(times),'.-');
    title('log10(Eval times (usec)) vs log10(evals per call)');
  end
end

% compute n random values uniformly distributed between min and max
function r = randrange(min, max, n)
  r = min + (max - min) * rand(n, 1);
end

% Validate certain parameters to avoid numeric and convergence problems
function tire = validate_parameters(tire, path)
  
  % Make sure we have reasonable values for certain fields
  if tire.VERTICAL_STIFFNESS == 0
    fprintf('%s: Invalid VERTICAL_STIFFNESS, using default\n', path);
    tire.VERTICAL_STIFFNESS = 200000;
  end
  if tire.Q_RE0 < .1
    fprintf('%s: Invalid Q_RE0, using default\n', path);
    tire.Q_RE0 = 1;
  end
  if tire.PDX1 < .01
    fprintf('%s: Invalid PDX1, using default\n', path);
    tire.PDX1 = 1;
  end
  if tire.PKX1 < 1
    fprintf('%s: Invalid PKX1, using default\n', path);
    tire.PKX1 = 25;
  end
  
  % MFeval 4.0 does not support tire side-flipping,
  % so we convert all tires to left tires here.
  if isprop(tire, 'TYRESIDE') && tire.TYRESIDE
    fprintf('MFeval does not support RIGHT tires: changing to LEFT\n');
  end
  tire.TYRESIDE = 0;
  
end

function compare_all()
  mftire_mfeval('tires/reference_tire_tvd3.tir');
  mftire_mfeval('tires/MagicFormula62_Parameters.tir');
  mftire_mfeval('tires/FSAE_Defaults.tir');
  mftire_mfeval('tires/PacejkaBook_Defaults.tir');
  mftire_mfeval('tires/MagicFormula52_Parameters.tir');
  mftire_mfeval('tires/MagicFormula61_Parameters.tir');
  
  % Model_Unstable.tir doesn't have longitudinal coefs: nothing we can do
  %%%mftire_mfeval('Model_Unstable.tir'); %% No longitudinal coefficients
end

% mftire_loadTIR.m - Load a user specified TIR file into a structure for use with MFeval.m
%
% Syntax: [CurrentTyre] = loadTIR(FileNameLocation)
% where FileNameLocation is a string including full path and name
%
% Based on MFeval 1.5 loadTIR.m
% Modified to add missing parameters, etc (see end of file)

% part of mftire 1.1.0
function originalTyre = mftire_loadTIR(FileNameLocation)
  
% NK: Based on mfeval 1.5 loadTIR.m version
% NK: modifications to add missing parameters, etc (see end of file)

fileID = fopen(FileNameLocation,'r'); % Open file filename for reading ('r')
% NK error handling
if fileID < 0
  fprintf('Can''t open %s\n', FileNameLocation);
  originalTyre = [];
  return;
end

counter = 0;
LineNumber = 0;

%While not at end of file
while feof(fileID) == 0
    
    currentLine = strtrim(fgetl(fileID));
    LineNumber = LineNumber + 1;
    
    if ~isempty(currentLine)
        %If current line is not a title, comment or description
        if (currentLine(1) ~= '[' && currentLine(1) ~= '$' && currentLine(1) ~= '!')
            %Append parameter counter
            counter = counter + 1;
            
            %Find location of = and $ in each row
            index1 = strfind(currentLine, '=') + 1;
            index2 = strfind(currentLine, '$') - 1;
            
            %Parameter Name
            value = strtrim(upper(currentLine(1:index1-2)));
            MFParamName{counter,1} = value;
            %Parameter Value
            if (isempty(index2) == 1)
                value = currentLine(index1:end);
                MFParamValue{counter,1} = strtrim(value);
            else
                value = currentLine(index1:index2);
                MFParamValue{counter,1} = strtrim(value);
            end
            %Parameter Description
            value = currentLine(index2+1:end);
            MFParamDescription{counter,1} = strtrim(value);
        end
    end
    
end

fclose(fileID);

clearvars -except MFParamName MFParamValue MFParamDescription FileNameLocation

% Mass is used as a variable name twice in the MF6.1 TIR
% format. Once in the Units section and once in the Inertia
% section. Rename the inertia section as MASS1. Undo this when
% writing the TIR file.
indSecondMass = find(strcmp(MFParamName,'MASS'),1,'last');
MFParamName(indSecondMass) = {'MASS1'};


%generate structure containing data of CurrentTyre
for ii = 1:length(MFParamName)
    if cellfun(@length,MFParamName(ii)) > 1
        if any(isstrprop(MFParamValue{ii},'digit')) > 0
            originalTyre.(MFParamName{ii}) = str2double(MFParamValue{ii});
            %fprintf('%16s = %.12g;\n', MFParamName{ii}, originalTyre.(MFParamName{ii}));
        else
            originalTyre.(MFParamName{ii}) = MFParamValue{ii};
            %fprintf('%16s = %s;\n', MFParamName{ii}, originalTyre.(MFParamName{ii}));
        end
    end
    
end
% include the description information in a separate field of
% the structure
for jj = 1:length(MFParamName)
    if cellfun(@length,MFParamName(jj)) > 1
        originalTyre.TIRDescription.(MFParamName{jj}) = MFParamDescription{jj};
    end
end

% NK: Following code added to ensure required parameters are present.

% Fields that might be missing and their default values
missing = {
  'NOMPRES'   0
  'INFLPRES'  0
  'PRESMIN'   0
  'PRESMAX'   0
  'PFZ1'      0
  'Q_RE0'     1
  'Q_V1'      0
  'Q_V2'      0
  'PPX1'      0
  'PPX2'      0
  'PPX3'      0
  'PPX4'      0
  'LKYC'      1
  'LKZC'      1
  'PEY5'      0
  'PKY4'      0
  'PKY5'      0
  'PKY6'      0
  'PKY7'      0
  'PPY1'      0
  'PPY2'      0
  'PPY3'      0
  'PPY4'      0
  'PPY5'      0
  'QDZ10'     0
  'QDZ11'     0
  'PPZ1'      0
  'PPZ2'      0
  'RBX3'      0
  'RBY4'      0
  'QSX4'      0
  'QSX5'      0
  'QSX6'      0
  'QSX7'      0
  'QSX8'      0
  'QSX9'      0
  'QSX10'     0
  'QSX11'     0
  'QSX12'     0
  'QSX13'     0
  'QSX14'     0
  'PPMX1'     0
  'QSY5'      0
  'QSY6'      0
  'QSY7'      0
  'QSY8'      0
  'Q_FZ2'     0
  'Q_FCX'     0
  'Q_FCY'     0
  'Q_RA1'     0
  'Q_RA2'     0
  'Q_RB1'     0
  'Q_RB2'     0
  'BOTTOM_OFFST' 0
  'Q_A1'      0
  'Q_A2'      0
  'PCFX1'     0
  'PCFX2'     0
  'PCFX3'     0
  'PCFY1'     0
  'PCFY2'     0
  'PCFY3'     0
  'FZMIN'     0
  'FZMAX'     0
  'ALPMIN'    -90*pi/180
  'ALPMAX'    90*pi/180
  'CAMMIN'    -8*pi/180
  'CAMMAX'    8*pi/180
  'KPUMIN'    -1
  'KPUMAX'    1
  'VERTICAL_STIFFNESS'      200000
  'LONGITUDINAL_STIFFNESS'  400000
  'LATERAL_STIFFNESS'       100000
  'YAW_STIFFNESS'           5000
  };

  for i=1:size(missing,1)
    name = missing{i,1};
    default = missing{i,2};
    if ~isfield(originalTyre, name)
      originalTyre.(name) = default;
    end
  end

  % Make sure we have reasonable values for certain fields
  if originalTyre.VERTICAL_STIFFNESS == 0
    fprintf('%s: Invalid VERTICAL_STIFFNESS, using default\n', FileNameLocation);
    originalTyre.VERTICAL_STIFFNESS = 200000;
  end
  
  % Make sure we have a reasonable load range
  if (originalTyre.FZMIN >= originalTyre.FZMAX)
    originalTyre.FZMIN = originalTyre.FNOMIN / 4;
    originalTyre.FZMAX = originalTyre.FNOMIN * 3;
  end

  if originalTyre.Q_RE0 < .1
    fprintf('%s: Invalid Q_RE0, using default\n', FileNameLocation);
    originalTyre.Q_RE0 = 1;
  end
  
  if originalTyre.PDX1 < .01
    fprintf('%s: Invalid PDX1, using default\n', FileNameLocation);
    originalTyre.PDX1 = 1;
  end

end
