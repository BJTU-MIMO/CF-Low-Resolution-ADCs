function [ alpha ] = finda( b )
switch b
    case{1}
        alpha=0.6366;
    case{2}
        alpha=0.8825;
    case{3}
        alpha=0.96546;
    case{4}
        alpha=0.990503;
    case{5}
        alpha=0.997501;
    otherwise
        alpha=1-(3.1415926*sqrt(3)/2)*2^(-2*b);
end
end

