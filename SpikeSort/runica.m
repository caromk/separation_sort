% runica() - Perform Independent Component Analysis (ICA) decomposition
%            of psychophysiological data using the ICA algorithm of 
%            Bell & Sejnowski (1995) with the natural gradient feature 
%            of Amari, Cichocki & Yang, or the extended-ICA algorithm 
%            of Lee, Girolami & Sejnowski (submitted).
% Usage:
%      simply >> [weights,sphere] = runica(data);
%       or    
%        else >> [weights,sphere,activations,bias,signs,lrates] ...
%                                 = runica(data,'Key1',Value1',...);
% Input_Variable:
%
% data        = input data (chans,frames*epochs). 
%               Note: each epoch should be baseline-zero'd using:
%                  >> data = rmbase(data,frames,basevector);
%
% Optional_Keywords       Keyword_Values                  Default_Values
%
% 'ncomps'    = number of ICA components to compute     (default -> chans)
% 'sphering'  = flag sphering of data ('on'/'off')      (default -> 'on')
% 'weights'   = initial weight matrix                   (default -> eye())
%                    (Note: if 'sphering' 'off',        (default -> spher())
% 'lrate'     = initial ICA learning rate << 1          (default -> heuristic)
% 'block'     = ICA block size (integer << datalength)  (default -> heuristic)
% 'anneal'    = annealing constant (0,1] (defaults -> 0.90, or 0.98, extended)
%                         controls speed of convergence
% 'annealdeg' = degrees weight change for annealing     (default -> 70)
% 'stop'      = stop training when weight-change < this (default -> 1e-6)
% 'maxsteps'  = maximum number of ICA training steps    (default -> 512)
% 'bias'      = perform bias adjustment ('on'/'off')    (default -> 'on')
% 'momentum'  = momentum gain [0,1]                     (default -> 0)
% 'extended'  = [N] perform tanh() "extended-ICA" with kurtosis estimation 
%               every N training blocks. If N < 0, fix number of sub-Gaussian
%               components to -N [faster than N>0] (default -> off)
% 'verbose'   = give ascii messages ('on'/'off')        (default -> 'on')
% 'order'     = order weights by variance ('on'/'off')  (default -> 'on')
%
% Output_Variables [RO = output in reverse order of projected mean variance 
%                        unless starting weight matrix passed ('weights' above)]
%
% weights     = ICA weight matrix (chans, comps)     [RO]
% sphere      = data sphering matrix (chans,chans) = spher(data)
%               Note: unmixing_matrix = weights*sphere {sphering off -> eye(chans)}
% activations = activation time courses of the output components (ncomps,frames*epochs)
% bias        = vector of final (ncomps) online bias [RO]    (default = zeros())
% signs       = extended-ICA signs for components    [RO]    (default = ones())
%                   [-1 = sub-Gaussian; 1 = super-Gaussian]
% lrates      = vector of learning rates used at each training step
%

% Toolbox Citation:
%
% Makeig, Scott. "ICA Toolbox for Psychophysiological Research (version 3.1.1)". 
% WWW Site, Computational Neurobiology Laboratory, The Salk Institute for Biological 
% Studies <www.cnl.salk.edu/~ica.html>, 1998. [World Wide Web Publication]. 
%
% Best Paper Reference:
%
% Scott Makeig, T-P. Jung, A. J. Bell, D. Ghahremani, T. J. Sejnowski,
% "Blind separation of auditory event-related brain responses into independent 
%  components," Proceedings of the National Academy of Sciences USA, 
% 94:10979-10984 (1997).
%
% First Reference:
%           Makeig, S., Bell, A.J., Jung, T-P and Sejnowski, T.J.,
%       "Independent component analysis of electroencephalographic data," 
%     In: D. Touretzky, M. Mozer and M. Hasselmo (Eds). Advances in Neural 
%  Information Processing Systems 8:145-151, MIT Press, Cambridge, MA (1996).
%
% For more information:
% http://www.cnl.salk.edu/~scott/icafaq.html - FAQ on ICA/EEG
% http://www.cnl.salk.edu/~scott/icabib.html - mss. on ICA & biosignals
% http://www.cnl.salk.edu/~tony/ica.html - math. mss. on ICA, with kernal code
%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Edit history %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  runica()  - by Scott Makeig, Tzyy-Ping Jung, Te-Won Lee and Tony Bell 
%                                  1996-97
%                       CNL/Salk Institute, La Jolla CA
%
%  04-30-96 from icatest.m and ~jung/.../wtwpwica.m -sm
%  07-28-97 new runica(), adds bias (default on), momentum (default off), 
%           extended-ICA (Lee & Sejnowski, 1997), cumulative angledelta 
%           (until lrate drops), keywords, signcount for speeding extended-ICA
%  10-07-97 put acos() outside verbose loop; verbose 'off' wasn't stopping -sm
%  11-11-97 adjusted help msg -sm
%  11-30-97 return eye(chans) if sphering 'off' or 'none' (undocumented option) -sm
%  02-27-98 use pinv() instead of inv() to rank order components when ncomps < chans -sm
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

function [weights,sphere,activations,bias,signs,lrates,y] = runica(data,p1,v1,p2,v2,p3,v3,p4,v4,p5,v5,p6,v6,p7,v7,p8,v8,p9,v9,p10,v10,p11,v11,p12,v12,p13,v13,p14,v14)

if nargin < 1
  help runica  
  return
end

[chans frames] = size(data); % determine the data size
datalength = frames;
if chans<2 | frames<chans
   fprintf('\nrunica() - data size (%d,%d) too small.\n\n', chans,frames);
   return
end
%
%%%%%%%%%%%%%%%%%%%%%% Declare defaults used below %%%%%%%%%%%%%%%%%%%%%%%%
%
MAX_WEIGHT           = 1e8;       % guess that weights larger than this have blown up
DEFAULT_STOP         = 0.000001;  % stop training if weight changes below this
DEFAULT_ANNEALDEG    = 60;        % when angle change reaches this value,
DEFAULT_ANNEALSTEP   = 0.90;      %     anneal by multiplying lrate by this
DEFAULT_EXTANNEAL    = 0.98;      %     or this if extended-ICA
DEFAULT_MAXSTEPS     = 512;       % ]top training after this many steps 
DEFAULT_MOMENTUM     = 0.0;       % default momentum weight

DEFAULT_BLOWUP       = 1000000000.0;   % = learning rate has 'blown up'
DEFAULT_BLOWUP_FAC   = 0.8;       % when lrate 'blows up,' anneal by this fac
DEFAULT_RESTART_FAC  = 0.9;       % if weights blowup, restart with lrate
                                  % lower by this factor
MIN_LRATE            = 0.000001;  % if weight blowups make lrate < this, quit
MAX_LRATE            = 0.1;       % guard against uselessly high learning rate
DEFAULT_LRATE        = 0.015/log(chans); 
                                  % heuristic default - may need adjustment
                                  %   for large or tiny data sets!
DEFAULT_BLOCK        = floor(sqrt(frames/3));  % heuristic default 
                                  % - may need adjustment!
% Extended-ICA option:
DEFAULT_EXTENDED     = 0;         % default off
DEFAULT_EXTBLOCKS    = 1;         % number of blocks per kurtosis calculation
DEFAULT_NSUB         = 1;         % initial default number of assumed sub-Gaussians
                                  % for extended-ICA
DEFAULT_EXTMOMENTUM  = 0.5;       % momentum term for computing extended-ICA kurtosis
MAX_KURTSIZE         = 6000;      % max points to use in kurtosis calculation
MIN_KURTSIZE         = 2000;      % minimum good kurtosis size (flag warning)
SIGNCOUNT_THRESHOLD  = 25;        % raise extblocks when sign vector unchanged
                                  % after this many steps
SIGNCOUNT_STEP       = 2;         % extblocks increment factor 

DEFAULT_SPHEREFLAG   = 'on';      % use the sphere matrix as the default
                                  %   starting weight matrix
DEFAULT_VERBOSE      = 1;         % write ascii info to calling screen
DEFAULT_BIASFLAG     = 1;         % default to using bias in the ICA update rule
DEFAULT_ORDER        = 1;         % default to ordering by variance 
%                                 
%%%%%%%%%%%%%%%%%%%%%%% Set up keyword default values %%%%%%%%%%%%%%%%%%%%%%%%%
%
if nargout < 2, 
    fprintf('runica() - needs at least two output arguments.\n');
    return
end
epochs = 1;							 % do not care how many epochs in data

sphering   = DEFAULT_SPHEREFLAG;     % default flags
verbose    = DEFAULT_VERBOSE;
order      = DEFAULT_ORDER;

block      = DEFAULT_BLOCK;          % heuristic default - may need adjustment!
lrate      = DEFAULT_LRATE;
annealdeg  = DEFAULT_ANNEALDEG;
annealstep = 0;                      % defaults declared below
nochange   = DEFAULT_STOP;
momentum   = DEFAULT_MOMENTUM;
maxsteps   = DEFAULT_MAXSTEPS;

weights    = 0;                      % defaults defined below
ncomps     = 0;
biasflag   = DEFAULT_BIASFLAG;

extended   = DEFAULT_EXTENDED;
extblocks  = DEFAULT_EXTBLOCKS;
kurtsize   = MAX_KURTSIZE;
signsbias  = 0.02;                   % bias towards super-Gaussian components
extmomentum= DEFAULT_EXTMOMENTUM;    % exp. average the kurtosis estimates
nsub       = DEFAULT_NSUB;
wts_blowup = 0;                      % flag =1 when weights too large
wts_passed = 0;                      % flag weights passed as argument
%
%%%%%%%%%% Collect keywords and values from argument list %%%%%%%%%%%%%%%
%
   if (nargin> 1 & rem(nargin,2) == 0)
      error('runica(): Even number of input arguments???')
   end
   for i = 3:2:nargin % for each Keyword
      Keyword = eval(['p',int2str((i-3)/2 +1)]);
      Value = eval(['v',int2str((i-3)/2 +1)]);
      if ~isstr(Keyword)
         error('runica(): keywords must be strings')
      end
      Keyword = lower(Keyword); % convert upper or mixed case to lower

      if strcmp(Keyword,'weights') | strcmp(Keyword,'weight')
         if isstr(Value)
            error(...
      'runica(): weights value must be a weight matrix or sphere')
         else
           weights = Value;
           wts_passed =1;
         end
      elseif strcmp(Keyword,'ncomps')
         if isstr(Value)
            error('runica(): ncomps value must be an integer')
         end
         ncomps = Value;
         if ~ncomps,
            ncomps = chans;
         end
      elseif strcmp(Keyword,'lrate')
         if isstr(Value)
            error('runica(): lrate value must be a number')
         end
         lrate = Value;
         if lrate>MAX_LRATE | lrate <0,
           error('runica(): lrate value is out of bounds'); 
         end
         if ~lrate,
            lrate = DEFAULT_LRATE;
         end
      elseif strcmp(Keyword,'block') | strcmp(Keyword,'blocksize')
         if isstr(Value)
            error('runica(): block size value must be a number')
         end
         block = Value;
         if ~block,
           block = DEFAULT_BLOCK; 
         end
      elseif strcmp(Keyword,'stop') | strcmp(Keyword,'nochange') ...
                    | strcmp(Keyword,'stopping')
         if isstr(Value)
            error('runica(): stop wchange value must be a number')
         end
         nochange = Value;
      elseif strcmp(Keyword,'maxsteps') | strcmp(Keyword,'steps')
         if isstr(Value)
            error('runica(): maxsteps value must be an integer')
         end
         maxsteps = Value;
         if ~maxsteps,
            maxsteps   = DEFAULT_MAXSTEPS;
         end
         if maxsteps < 0
            error('runica(): maxsteps value must be a positive integer')
         end
      elseif strcmp(Keyword,'anneal') | strcmp(Keyword,'annealstep')
         if isstr(Value)
            error('runica(): anneal step constant must be a number (0,1)')
         end
         annealstep = Value;
         if annealstep <=0 | annealstep > 1,
            error('runica(): anneal step value must be (0,1]')
         end
      elseif strcmp(Keyword,'annealdeg') | strcmp(Keyword,'degrees')
         if isstr(Value)
            error('runica(): annealdeg value must be a number')
         end
         annealdeg = Value;
         if ~annealdeg,
             annealdeg = DEFAULT_ANNEALDEG;
         elseif annealdeg > 180 | annealdeg < 0
          error('runica(): annealdeg value is out of bounds [0,180]')
                                              
         end
      elseif strcmp(Keyword,'momentum')
         if isstr(Value)
            error('runica(): momentum value must be a number')
         end
         momentum = Value;
         if momentum > 1.0 | momentum < 0
          error('runica(): momentum value is out of bounds [0,1]')
         end
      elseif strcmp(Keyword,'sphering') | strcmp(Keyword,'sphereing') ...
                | strcmp(Keyword,'sphere')
         if ~isstr(Value)
           error('runica(): sphering value must be on, off, or none')
         else 
           Value = lower(Value);
           if ~strcmp(Value,'on') & ~strcmp(Value,'off') & ~strcmp(Value,'none'),
             error('runica(): sphering value must be on or off')
           end
           sphering = Value;
         end
      elseif strcmp(Keyword,'bias')
         if ~isstr(Value)
           error('runica(): bias value must be on or off')
         else 
           Value = lower(Value);
           if strcmp(Value,'on') 
              biasflag = 1;
           elseif strcmp(Value,'off'),
              biasflag = 0;
           else
              error('runica(): bias value must be on or off')
           end
         end
      elseif strcmp(Keyword,'extended') | strcmp(Keyword,'extend')
         if isstr(Value)
           error('runica(): extended value must be an integer >=0')
         else
           extended = 1;      % turn on extended-ICA
           extblocks = fix(Value); % number of blocks per kurt() compute
           if extblocks < 0
                nsub = -1*fix(extblocks);  % fix this many sub-Gauss comps
           elseif ~extblocks,
                extended = 0;             % turn extended-ICA off
           elseif kurtsize>frames,   % length of kurtosis calculation
                kurtsize = frames;
                if kurtsize < MIN_KURTSIZE
                   fprintf(...
   'runica() warning: kurtosis values inexact for << %d points.\n',...
                                                         MIN_KURTSIZE);
                end
           end
         end
      elseif strcmp(Keyword,'verbose') 
         if ~isstr(Value)
            error('runica(): verbose flag value must be on or off')
         elseif strcmp(Value,'on'),
             verbose = 1; 
         elseif strcmp(Value,'off'),
             verbose = 0; 
         else
             error('runica(): verbose flag value must be on or off')
         end
      elseif strcmp(Keyword,'order') 
         if ~isstr(Value)
            error('runica(): order flag value must be on or off')
         elseif strcmp(Value,'on'),
             order = 1; 
         elseif strcmp(Value,'off'),
             order = 0; 
         else
             error('runica(): order flag value must be on or off')
         end
      else
         error('runica(): unknown flag')
      end
   end
%
%%%%%%%%%%%%%%%%%%%%%%%% Initialize weights, etc. %%%%%%%%%%%%%%%%%%%%%%%%
%
if ~annealstep,
  if ~extended,
    annealstep = DEFAULT_ANNEALSTEP;     % defaults defined above
  else
    annealstep = DEFAULT_EXTANNEAL;       % defaults defined above
  end
end % else use annealstep from commandline

if ~annealdeg, 
    annealdeg  = DEFAULT_ANNEALDEG - momentum*90; % heuristic
    if annealdeg < 0,
        annealdeg = 0;
    end
end
if ~ncomps,
    ncomps = chans;
elseif ncomps >  chans | ncomps < 1
    fprintf('runica(): number of components must be 1 to %d.\n',chans);
    return
end

if weights ~= 0,                    % initialize weights
  % starting weights are being passed to runica() from the commandline
    if verbose,
       fprintf('Using starting weight matrix named in argument list ...\n')
    end
    if  chans>ncomps & weights ~=0,
        [r,c]=size(weights);
        if r~=ncomps | c~=chans,
     fprintf(...
      'runica(): weight matrix must have %d rows, %d columns.\n', ...
                              chans,ncomps);
            return;
        end
    end
end;   
%
%%%%%%%%%%%%%%%%%%%%% Check keyword values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if frames<chans,
    fprintf('runica(): data length %d < data channels %f!\n',frames,chans)
    return
elseif block < 2,
    fprintf('runica(): block size %d too small!\n',block)
    return
elseif block > frames, 
    fprintf('runica(): block size exceeds data length!\n');
    return
elseif floor(epochs) ~= epochs,
    fprintf('runica(): data length is not a multiple of the epoch length!\n');
    return
elseif nsub > ncomps
    fprintf('runica(): there can be at most %d sub-Gaussian components!\n',ncomps);
    return
end;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Process the data %%%%%%%%%%%%%%%%%%%%%%%%%%
%
if verbose,
   fprintf( ...
            '\nInput data size [%d,%d] = %d channels, %d frames.\n', ...
                        chans,frames,chans,frames);
   if ~extended
      fprintf('Finding %d ICA components using logistic ICA.\n',ncomps);
   else % if extended
      fprintf('Finding %d ICA components using extended ICA.\n',ncomps);
      if extblocks > 0
         fprintf(...
 'Kurtosis will be calculated initially every %d blocks using %d data points.\n',...
                                            extblocks,     kurtsize);
      else
         fprintf(...
 'Kurtosis will not be calculated. Exactly %d sub-Gaussian components assumed.\n',...
                                        nsub);
      end
   end
   fprintf('Initial learning rate will be %g, block size %d.\n',lrate,block);
   if momentum>0,
      fprintf('Momentum will be %g.\n',momentum);
   end
   fprintf( ...
'Learning rate will be multiplied by %g whenever angledelta >= %g deg.\n', ...
               annealstep,annealdeg);
 fprintf('Training will end when wchange < %g or after %d steps.\n', ...
                    nochange,maxsteps);
  if biasflag,
    fprintf('Online bias adjustment will be used.\n');
  else
    fprintf('Online bias adjustment will not be used.\n');
  end
end
%
%%%%%%%%%%%%%%%%%%%%%%%%% Remove overall row means %%%%%%%%%%%%%%%%%%%%%%%%
%
if verbose,
    fprintf('Removing mean of each channel ...\n');
end
data = data - mean(data')'*ones(1,frames);      % subtract row means

if verbose,
   fprintf('Final training data range: %g to %g\n', ...
                          min(min(data)),max(max(data)));
end
%
%%%%%%%%%%%%%%%%%%% Use the sphering matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

if strcmp(sphering,'on'), %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if verbose,
      fprintf('Computing the sphering matrix...\n');
  end
  sphere = 2.0*inv(sqrtm(cov(data'))); % find the "sphering" matrix = spher()
  if ~weights,
      if verbose,
          fprintf('Starting weights are the identity matrix ...\n');
      end
      weights = eye(ncomps,chans); % begin with the identity matrix
  else % weights given on commandline
      if verbose,
          fprintf('Using starting weights named on commandline ...\n');
      end
  end
  if verbose,
      fprintf('Sphering the data ...\n');
  end
  data = sphere*data;      % actually decorrelate the electrode signals

elseif strcmp(sphering,'off') %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if ~weights
      if verbose,
       fprintf('Using the sphering matrix as the starting weight matrix ...\n');
       fprintf('Returning the identity matrix in variable "sphere" ...\n');
      end
      sphere = 2.0*inv(sqrtm(cov(data'))); % find the "sphering" matrix = spher()
      weights = eye(ncomps,chans)*sphere;  % begin with the sphering matrix
      sphere = eye(chans);                 % return the identity matrix
  else % weights ~= 0
      if verbose,
       fprintf('Using starting weights named on commandline ...\n');
       fprintf('Returning the identity matrix in variable "sphere" ...\n');
      end
      sphere = eye(chans);                 % return the identity matrix
  end
elseif strcmp(sphering,'none')
  sphere = eye(chans);                     % return the identity matrix
  if ~weights
      if verbose,
       fprintf('Starting weights are the identity matrix ...\n');
       fprintf('Returning the identity matrix in variable "sphere" ...\n');
      end
      weights = eye(ncomps,chans); % begin with the identity matrix
  else % weights ~= 0
      if verbose,
       fprintf('Using starting weights named on commandline ...\n');
       fprintf('Returning the identity matrix in variable "sphere" ...\n');
      end
  end
  sphere = eye(chans,chans);
  if verbose,
      fprintf('Returned variable "sphere" will be the identity matrix.\n');
  end
end
%
%%%%%%%%%%%%%%%%%%%%%%%% Initialize ICA training %%%%%%%%%%%%%%%%%%%%%%%%%
%
  lastt=fix((datalength/block-1)*block+1);
  BI=block*eye(ncomps,ncomps);
  delta=zeros(1,chans*ncomps);
  changes = [];
  degconst = 180./pi;
  startweights = weights;
  prevweights = startweights;
  oldweights = startweights;
  prevwtchange = zeros(chans,ncomps);
  oldwtchange = zeros(chans,ncomps);
  lrates = zeros(1,maxsteps);
  onesrow = ones(1,block);
  bias = zeros(ncomps,1);
  signs = ones(1,ncomps);    % initialize signs to nsub -1, rest +1
  for k=1:nsub
      signs(k) = -1;
  end
  if extended & extblocks < 0 & verbose,
    fprintf('Fixed extended-ICA sign assignments:  ');
    for k=1:ncomps
       fprintf('%d ',signs(k));	
    end; fprintf('\n');
  end
  signs = diag(signs); % make a diagonal matrix
  oldsigns = zeros(size(signs));;
  signcount = 0;              % counter for same-signs
  signcounts = [];
  urextblocks = extblocks;    % original value, for resets
  old_kk = zeros(1,ncomps);   % for kurtosis momemtum
%
%%%%%%%% ICA training loop using the logistic sigmoid %%%%%%%%%%%%%%%%%%%
%
  if verbose,
    fprintf('Beginning ICA training ...');
    if extended,
        fprintf(' first training step may be slow ...\n');
    else
        fprintf('\n');
    end
  end
  step=0;
  laststep=0; 
  blockno = 1;  % running block counter for kurtosis interrupts

  while step < maxsteps, %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    permute=randperm(datalength); % shuffle data order at each step

    for t=1:block:lastt, %%%%%%%%% ICA Training Block %%%%%%%%%%%%%%%%%%%
      if biasflag                                                   
       u=weights*data(:,permute(t:t+block-1)) + bias*onesrow;      
      else                                                             
       u=weights*data(:,permute(t:t+block-1));                      
      end                                                              
      if ~extended
       %%%%%%%%%%%%%%%%%%% Logistic ICA weight update %%%%%%%%%%%%%%%%%%%
       y=1./(1+exp(-u));                                                %
       weights=weights+lrate*(BI+(1-2*y)*u')*weights;                   %
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      else % extended-ICA
       %%%%%%%%%%%%%%%%%%% Extended-ICA weight update %%%%%%%%%%%%%%%%%%%
       y=tanh(u);                                                       %
       weights = weights + lrate*(BI-signs*y*u'-u*u')*weights;          %
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      end
      if biasflag 
         if ~extended
       %%%%%%%%%%%%%%%%%%%%%%%% Logistic ICA bias %%%%%%%%%%%%%%%%%%%%%%%
             bias = bias + lrate*sum((1-2*y)')'; % for logistic nonlin. %
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         else % extended
       %%%%%%%%%%%%%%%%%%% Extended-ICA bias %%%%%%%%%%%%%%%%%%%%%%%%%%%%
             bias = bias + lrate*sum((-2*y)')';  % for tanh() nonlin.   %
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         end                                    
      end

      if momentum > 0 %%%%%%%%% Add momentum %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        weights = weights + momentum*prevwtchange;                
        prevwtchange = weights-prevweights;                      
        prevweights = weights;                                  
      end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      if max(max(abs(weights))) > MAX_WEIGHT
          wts_blowup = 1;
          change = nochange;
      end
      if extended & ~wts_blowup
        %%%%%%%%%% Extended-ICA kurtosis estimation %%%%%%%%%%%%%%%%%%%%%
       if extblocks > 0 & rem(blockno,extblocks) == 0, 
                                  % recompute signs vector using kurtosis
        if kurtsize < frames 
            rp = randperm(frames);             	    % random data subset
            partact=weights*data(:,rp(1:kurtsize)); % partial activation
        else                                        % for small data sets,
            partact=weights*data;                   % use whole data
        end
        m2=mean(partact'.^2).^2;             
        m4= mean(partact'.^4);
        kk= (m4./m2)-3.0;                           % kurtosis estimates
        if extmomentum
            kk = extmomentum*old_kk + (1.0-extmomentum)*kk; % use momentum
            old_kk = kk;
        end
        signs=diag(sign(kk+signsbias));             % pick component signs
        if signs == oldsigns,
           signcount = signcount+1;
        else
           signcount = 0;
        end
        oldsigns = signs;
        signcounts = [signcounts signcount];
        if signcount >= SIGNCOUNT_THRESHOLD,
         extblocks = fix(extblocks * SIGNCOUNT_STEP);% make kurt() estimation
         signcount = 0;                             % less frequent if sign
        end                                         % is not changing
       end % extblocks > 0 & . . .
      end % if extended %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      blockno = blockno + 1;
      if wts_blowup
         break
      end
    end % training block %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~wts_blowup
      oldwtchange = weights-oldweights;
      step=step+1; 
      %
      %%%%%%% Compute and print weight and update angle changes %%%%%%%%%
      %
      lrates(1,step) = lrate;
      angledelta=0.;
      delta=reshape(oldwtchange,1,chans*ncomps);
      change=delta*delta'; 
    end
    %
    %%%%%%%%%%%%%%%%%%%%%% Restart if weights blow up %%%%%%%%%%%%%%%%%%%%
    %
    if wts_blowup | isnan(change)|isinf(change),  % if weights blow up,
        %fprintf('');
        step = 0;                          % start again
        change = nochange;
        wts_blowup = 0;                    % re-initialize variables
        blockno = 1;
        lrate = lrate*DEFAULT_RESTART_FAC; % with lower learning rate
        weights = startweights;            % and original weight matrix
        oldweights = startweights;            
        change = nochange;
        oldwtchange = zeros(chans,ncomps);
        delta=zeros(1,chans*ncomps);
        olddelta = delta;
        extblocks = urextblocks;
        prevweights = startweights;
        prevwtchange = zeros(chans,ncomps);
        lrates = zeros(1,maxsteps);
        bias = zeros(ncomps,1);
        if extended
           signs = ones(1,ncomps);    % initialize signs to nsub -1, rest +1
           for k=1:nsub
               signs(k) = -1;
           end
           signs = diag(signs); % make a diagonal matrix
           oldsigns = zeros(size(signs));;
        end
        if lrate> MIN_LRATE
          r = rank(data);
          if r<ncomps
           fprintf('Data has rank %d. Cannot compute %d components.\n',...
                                  r,ncomps);
           return
          else
           if verbose,
            fprintf(...
             'Lowering learning rate to %g and starting again.\n',lrate);
           end
          end
        else
          fprintf( ...
          'runica(): QUITTING - weight matrix may not be invertible!\n');
          return;
        end
    else % if weights in bounds 
     %
     %%%%%%%%%%%%% Print weight update information %%%%%%%%%%%%%%%%%%%%%%
     %
     if step> 2 
         angledelta=acos((delta*olddelta')/sqrt(change*oldchange));
     end
     if verbose,
      if step > 2, 
          if ~extended,
              fprintf(...
         'step %d - lrate %5f, wchange %7.6f, angledelta %4.1f deg\n', ...
          step,lrate,change,degconst*angledelta);
          else
              fprintf(...
'step %d - lrate %5f, wchange %7.6f, angledelta %4.1f deg, %d subgauss\n',...
          step,lrate,change,degconst*angledelta,(ncomps-sum(diag(signs)))/2);
          end
      elseif ~extended
         fprintf(...
          'step %d - lrate %5f, wchange %7.6f\n',step,lrate,change);
      else
         fprintf(...
          'step %d - lrate %5f, wchange %7.6f, %d subgauss\n',...
           step,lrate,change,(ncomps-sum(diag(signs)))/2);
      end % step > 2
     end; % if verbose
%
%%%%%%%%%%%%%%%%%%%% Save current values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
      changes = [changes change];
      oldweights = weights;
%
%%%%%%%%%%%%%%%%%%%% Anneal learning rate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
      if degconst*angledelta > annealdeg,  
        lrate = lrate*annealstep;          % anneal learning rate
        olddelta   = delta;                % accumulate angledelta until
        oldchange  = change;               %  annealdeg is reached
      elseif step == 1                     % on first step only
        olddelta   = delta;                % initialize 
        oldchange  = change;               
      end
%
%%%%%%%%%%%%%%%%%%%% Apply stopping rule %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
      if step >2 & change < nochange,      % apply stopping rule
        laststep=step;            
           step=maxsteps;                  % stop when weights stabilize
      elseif change > DEFAULT_BLOWUP,      % if weights blow up,
        lrate=lrate*DEFAULT_BLOWUP_FAC;    % keep trying 
      end;                                 % with a smaller learning rate
    end; % end if weights in bounds

  end; % end training %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if ~laststep
    laststep = step;
  end;
  lrates = lrates(1,1:laststep);           % truncate lrate history vector
%
%%%%%% Sort components in descending order of max projected variance %%%%
%
if verbose,
  fprintf(...
 'Sorting components in descending order of mean projected variance ...\n');
end

if wts_passed == 0 & order == 1
  activations=weights*data; % filter the data -> activation wave forms
  %
  %%%%%%%%%%%%%%%%%%%% Find mean variances %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %

  meanvar  = zeros(ncomps,1);      % size of the projections
  if ncomps == chans % if weights are square . . .
    winv = inv(weights*sphere);
  else
    fprintf('Using pseudo-inverse of weight matrix to rank order component projections.\n');
    winv = pinv(weights*sphere);
  end
  for s=1:ncomps
    if verbose,
      fprintf('%d ',s);         % construct single-component data matrix
    end
   					% project to scalp, then add row means 
    compproj = winv(:,s)*activations(s,:);
    meanvar(s) = mean(sum(compproj.*compproj)/(size(compproj,1)-1));
                                            % compute mean variance 
  end                                         % at all scalp channels
  if verbose,
   fprintf('\n');
  end
  %
  %%%%%%%%%%%%%% Sort components by mean variance %%%%%%%%%%%%%%%%%%%%%%%%
  %
  [sortvar, windex] = sort(meanvar);
  windex = windex(ncomps:-1:1); % order large to small 
  meanvar = meanvar(windex);
  % 
  %%%%%%%%%%%%%%%%%%%%% Filter data using final weights %%%%%%%%%%%%%%%%%%
  %
  if nargout>2, % if activations are to be returned
   if verbose,
     fprintf('Permuting the activation wave forms ...\n');
   end
   activations = activations(windex,:);
  else
   clear activations
  end
  weights = weights(windex,:);% reorder the weight matrix
  bias  = bias(windex);		% reorder them
  signs = diag(signs);        % vectorize the signs matrix
  signs = signs(windex);      % reorder them
else
  if verbose,
   fprintf('Components not ordered by variance.\n');
  end
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
    
return

if nargout > 6
  u=weights*data + bias*ones(1,frames);      
  y = zeros(size(u));
  for c=1:chans
    for f=1:frames
       y(c,f) = 1/(1+exp(-u(c,f)));
    end
  end
end
