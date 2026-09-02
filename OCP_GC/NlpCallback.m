%
%     This file is part of CasADi.
%
%     CasADi -- A symbolic framework for dynamic optimization.
%     Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
%                             K.U. Leuven. All rights reserved.
%     Copyright (C) 2011-2014 Greg Horn
%
%     CasADi is free software; you can redistribute it and/or
%     modify it under the terms of the GNU Lesser General Public
%     License as published by the Free Software Foundation; either
%     version 3 of the License, or (at your option) any later version.
%
%     CasADi is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%     Lesser General Public License for more details.
%
%     You should have received a copy of the GNU Lesser General Public
%     License along with CasADi; if not, write to the Free Software
%     Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
%
%

classdef NlpCallback < casadi.Callback
  properties
    data
    nx
    ng
    iter_count = 0;
  end
  methods
    function self = NlpCallback(name, nx, ng)
      self@casadi.Callback();
      self.nx = nx;
      self.ng = ng;
      construct(self, name);
      self.data = [];
    end
    function [returncode] = eval(self, arg)
        x_val = full(arg{1});
        
        % Incremental counter
        % if ~isfield(self, 'iter_count')
        %     self.iter_count = 1;
        % else
            self.iter_count = self.iter_count + 1;
        % end

        self.data = [ self.data full(arg{1})];
      
        % Save the first 20 directly
        if self.iter_count <= 20
            filename = sprintf('first_iter_%04d.mat', self.iter_count);
            x_snapshot = x_val; %#ok<NASGU>
            save(filename, 'x_snapshot');
        end
        if any(self.iter_count == [250:270 620:656])
            filename = sprintf('first_iter_%04d.mat', self.iter_count);
            x_snapshot = x_val; %#ok<NASGU>
            save(filename, 'x_snapshot');
        end
        % % Save a rolling log of last 20 iterations (overwrite files)
        % elseif self.iter_count > 2890
        %     filename = sprintf('last_iter_%04d.mat', self.iter_count);
        %     x_snapshot = x_val;
        %     save(filename, 'x_snapshot');
        % end

      returncode = {0};
    end
    
    function out = get_sparsity_in(self,i)
      n = casadi.nlpsol_out(i);
      if strcmp(n,'f')
        out = [1 1];
      elseif strcmp(n,'lam_x') || strcmp(n,'x')
        out = [self.nx 1];
      elseif strcmp(n,'lam_g') || strcmp(n,'g')
        out = [self.ng 1];
      else
        out = [0 0];
      end
      out = casadi.Sparsity.dense(out(1),out(2));
    end
    function out = get_n_in(self)
      out = casadi.nlpsol_n_out();
    end
    function out = get_n_out(self)
      out = 1;
    end
  end
end