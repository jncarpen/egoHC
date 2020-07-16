function G = gauss2D(x, y, sigma)
%2DGAUSSIAN Make 2-dimensional gaussian filter
%   fix this function up to be more useful *

%   USAGE:
%   G = gauss2D(-5:5, -5:5, 1);

for xx = 1:length(x)
    for yy = 1:length(y)
        G(xx,yy) = (1/(2*pi*sigma^2))*exp(-((x(xx)^2)+(y(yy)^2))/(2*sigma^2));
    end
end

end

