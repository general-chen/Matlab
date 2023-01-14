function [Fx_smoothed] = smoother200(Fx)
    
   Fx_smoothed=smooth(Fx,200,'sgolay'); 

end

