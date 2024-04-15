function G = LinearPathProfile(pMat, l1, l2, Tmax)
% pMat: matrix of points: first col: x coordinates, second col: y
% coordinates : each row represents a point of the form [x, y].
% a point is always added to the last row of the matrix to close the loop.


nPoints = size(pMat, 1);
nPaths = nPoints;

intervals = Tmax/nPaths * (0:1:nPaths);
timeSubIntervals = Tmax/nPaths;

% close the loop
pMat = [pMat;pMat(1,:)];

[f, fp] = UnitRise7();


    function [phi1, phi1p, phi2, phi2p] = F(t)
        phi1 = zeros(1, numel(t));
        phi1p = zeros(1, numel(t));
        phi2 = zeros(1, numel(t));
        phi2p = zeros(1, numel(t));

        N = numel(t);
        for k=1:N

            if (t(k)<0 || t(k)>=Tmax)

                A = pMat(end-1, :);
                B = pMat(end,:);
                nhat = (B-A);

                lambda = f(timeSubIntervals/timeSubIntervals);
                lambdap = 1/timeSubIntervals * fp(timeSubIntervals/timeSubIntervals);

                X = A + lambda*nhat;
                C1 = (l1^2 + l2^2);
                C2 = (2*l1*l2);

                X1p = nhat(1)*lambdap;
                X2p = nhat(2)*lambdap;

                D = (vecnorm(X)^2 - C1)/C2; % = cos(phi2)
                phi2(k) = -acos(D);
                phi2p(k) = 1/sqrt(1-D^2) * (1/C2 * 2*(X(1)*X1p + X(2)*X2p));

                phi1(k) = atan2(X(2), X(1)) - atan2(l2*sin(phi2(k)), l1+l2*D);

                den = (l2*sin(phi2(k)))^2 + (l1+l2*D)^2;
                der1 = -l2 * sin(phi2(k)) * phi2p(k);
                der2 = l2 * cos(phi2(k)) * phi2p(k);

                phi1p(k) = (-X(2)/(vecnorm(X)^2) * X1p) + (X(1)/(vecnorm(X)^2) * X2p) + ...
                    -(...
                    (-(l2*sin(phi2(k)))/den * der1) + ((l1+l2*D)/den * der2) ...
                    );

            else
                for i=1:nPoints
                    if(t(k)>=intervals(i) && t(k)<intervals(i+1))

                        A = pMat(i, :);
                        B = pMat(i+1,:);
                        nhat = (B-A);


                        lambda = f(mod(t(k),timeSubIntervals)/timeSubIntervals);
                        lambdap = 1/timeSubIntervals * fp(mod(t(k),timeSubIntervals)/timeSubIntervals);

                        X = A + lambda*nhat;
                        C1 = (l1^2 + l2^2);
                        C2 = (2*l1*l2);

                        X1p = nhat(1)*lambdap;
                        X2p = nhat(2)*lambdap;

                        D = (vecnorm(X)^2 - C1)/C2; % = cos(phi2)
                        phi2(k) = -acos(D);
                        phi2p(k) = 1/sqrt(1-D^2) * (1/C2 * 2*(X(1)*X1p + X(2)*X2p));

                        phi1(k) = atan2(X(2), X(1)) - atan2(l2*sin(phi2(k)), l1+l2*D);

                        den = (l2*sin(phi2(k)))^2 + (l1+l2*D)^2;
                        der1 = -l2 * sin(phi2(k)) * phi2p(k);
                        der2 = l2 * cos(phi2(k)) * phi2p(k);

                        phi1p(k) = (-X(2)/(vecnorm(X)^2) * X1p) + (X(1)/(vecnorm(X)^2) * X2p) + ...
                            -(...
                            (-(l2*sin(phi2(k)))/den * der1) + ((l1+l2*D)/den * der2) ...
                            );
                        break;
                    end
                end
            end
        end
    end

G = @F;
end

%
% lambda = f(t(k)/Tmax);
%                     lambdap = 1/Tmax * fp(t(k)/Tmax);
%
%                     X = A + lambda*nhat;
%                     C1 = (l1^2 + l2^2);
%                     C2 = (2*l1*l2);
%
%                     X1p = nhat(1)*lambdap;
%                     X2p = nhat(2)*lambdap;
%
%                     D = (vecnorm(X)^2 - C1)/C2; % = cos(phi2)
%                     phi2(k) = -acos(D);
%                     phi2p(k) = 1/sqrt(1-D^2) * (1/C2 * 2*(X(1)*X1p + X(2)*X2p));
%
%                     phi1(k) = atan2(X(2), X(1)) - atan2(l2*sin(phi2(k)), l1+l2*D);
%
%                     den = (l2*sin(phi2(k)))^2 + (l1+l2*D)^2;
%                     der1 = -l2 * sin(phi2(k)) * phi2p(k);
%                     der2 = l2 * cos(phi2(k)) * phi2p(k);
%
%                     phi1p(k) = (-X(2)/(vecnorm(X)^2) * X1p) + (X(1)/(vecnorm(X)^2) * X2p) + ...
%                         -(...
%                         (-(l2*sin(phi2(k)))/den * der1) + ((l1+l2*D)/den * der2) ...
%                         );
