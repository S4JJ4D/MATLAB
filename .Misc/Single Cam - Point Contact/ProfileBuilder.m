function G = ProfileBuilder(p)
% p is an array of the form:
% p = [t1 x1 t2 x2 ...]

if mod(length(p),2) ~= 0
    error('length(p) must be an even number');
end

nPoints = length(p)/2;


[f, fp, fpp] = UnitRise7();

    function [x, xdot, xddot] = F(t)
        x = zeros(1, numel(t));
        xdot = zeros(1, numel(t));
        xddot = zeros(1, numel(t));
        N = numel(t);
        for k=1:N
            if (t(k) >= p(end-1))
                x(k) = p(end);
                xdot(k) = 0;
                xddot(k) = 0;
            elseif (t(k) < p(1))
                x(k) = p(2);
                xdot(k) = 0;
                xddot(k) = 0;
            else
                for i=1:nPoints
                    if(t(k)>=p(2*i-1) && t(k)<p(2*i+1))
                        delta_t = p(2*i+1) - p(2*i-1);
                        delta_x = p(2*i+2) - p(2*i);
                        x(k) = delta_x * f( (t(k)-p(2*i-1))/delta_t ) + p(2*i);
                        xdot(k) = (delta_x/delta_t) *fp( (t(k)-p(2*i-1))/delta_t );
                        xddot(k) = (delta_x/(delta_t.^2)) *fpp( (t(k)-p(2*i-1))/delta_t );
                        break;
                    end
                end
            end
        end
    end


G = @F;

end