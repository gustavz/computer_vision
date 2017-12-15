function [Fx,Fy] = sobel_xy(Image)
% In dieser Funktion soll das Sobel-Filter implementiert werden, welches
% ein Graustufenbild einliest und den Bildgradienten in x- sowie in
% y-Richtung zurückgibt.

% Interpolationsfilter
di = [1  2 1];

% Ableitungsfilter
dd = [1 0 -1];

% Ausnutzen der Separabilität des Sobel-Filters. 
Fx=conv2(di,dd,Image,'same');
Fy=conv2(dd,di,Image,'same');

end

