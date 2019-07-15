classdef Opt
    
    properties
        opt
    end
    
    methods
        function O = Opt(vargin) % vargin is intentionally not varargin - takes a structure
            n = length(vargin);
            opt = struct();
            if (ceil(n/2) ~= n/2) % check if even
                error('Opt:Opt','Opt requires an even number of arguments');
            else
                for i = 1:2:n
                    O.opt = setfield(O.opt,vargin{i},vargin{i+1});
                end
            end
        end
        
        function value = get(self,key,varargin)
            if isfield(self.opt,key)
                value = getfield(self.opt,key);
            elseif ~isempty(varargin)
                value = varargin{1};
            else
                value = [];
            end
        end
    end
                    
end