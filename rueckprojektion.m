function repro_error=rueckprojektion(Korrespondenzen, P1, Image2, T, R, K)

% Anzahl der KP-Paare
N = size(Korrespondenzen,2);

% Berechne die Weltkoordinaten der 3D-Punkte (in Bild1 berechnet) in Kamera 2
P1_2 = bsxfun(@plus, R*P1, T);

% Berechne die Rückprojektion der 3D-Punkte in Bild 2
x2_reprojected = K * bsxfun(@times, P1_2, 1./P1_2(3,:));

% Stelle die Korrespondenzpunktpaare dar und zeichne die Rückprojektion der
% 3D-Punkte in Bild 2 ein
figure('name', 'Rückprojektion von Korrespondenzpunkten');
imshow(Image2, 'Border','tight');

% Initialisierung des Rückprojektionsfehlers
repro_error=0;

for i = 1:N
    hold on;
    % Punkte im zweiten Bild
    plot(Korrespondenzen(3,i), Korrespondenzen(4,i), 'Marker', 'x', 'Color', [0 1 0], 'MarkerSize', 3);
    text(Korrespondenzen(3,i), Korrespondenzen(4,i), sprintf('%d', i), 'fontsize', 10, 'color', [0 1 0]); 
    % Rückprojezierte 3-D-Punkte aus Bild1 in Bild2
    plot(x2_reprojected(1,i), x2_reprojected(2,i), 'Marker', 'x', 'Color', [1 0 0], 'MarkerSize',3);
    text(x2_reprojected(1,i), x2_reprojected(2,i), sprintf('%d', i), 'fontsize', 10, 'color', [1 0 0]);
    
    repro_error=repro_error + norm(x2_reprojected(1:2,i)-Korrespondenzen(3:4,i));
end
% Mittele den Rückprojektionsfehler über alle Punkte
repro_error=repro_error/N;
hold off;

end