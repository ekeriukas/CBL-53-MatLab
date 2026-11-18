function [V] = CylinderVolume(Ca,Cyl)
% This function provides the cylinder volume as function of 
% Ca : Crankangle [degrees]
% Cyl :  a struct containing
%   Cyl.S : Stroke
%   Cyl.B                   : Bore
%   Cyl.ConRod              : Connecting Rod length
%   Cyl.CompressionRatio    : Compession Ratio
%   Cyl.TDCangle            : Angle associated with the Top Dead Center
%----------------------------------------------------------------------
B   = Cyl.Bore;
S   = Cyl.Stroke;
cr  = Cyl.CompressionRatio;
r   = S/2;
l   = Cyl.ConRod;
%-------------------------------------------------------------------------------------------------------
CAl     = Ca-Cyl.TDCangle;
Vd      = pi*(B/2)^2*S;
Vc      = Vd/(cr-1);
Ap      = pi * (B/2)^2;
% Cylinder volume as a function of crank angle (CA in degrees)
V   = Vc + Ap .* (r .* (1 - cosd(CAl)) + l .* (1 - sqrt(1 - (r ./ l).^2 .* sind(CAl).^2)));





