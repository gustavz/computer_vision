function [T,R, lambda, P1] = rekonstruktion(T1,T2,R1,R2, Korrespondenzen, K)

%% Finden der korrekten Euklidischen Bewegung über die Positive-Tiefen-Bedingung


% Speichere die Ts und Rs in einem Array ab
T_all = cat(2, T1, T2, T1, T2);
R_all = cat(3, R1, R1, R2, R2);

N = length(Korrespondenzen);
x1_pixel_array = [Korrespondenzen(1:2,:); ones(1, N)];
x2_pixel_array = [Korrespondenzen(3:4,:); ones(1, N)];

% Kalibrierung der Pixelkoordinaten um Bildkoordinaten zu erhalten
x1_array = K \ x1_pixel_array;
x2_array = K \ x2_pixel_array;

% 4 mögliche Kombinationen aus T_1,2 und R_1,2, d.h. wir muessen 4 mal das
% Gleichungssystem M1 aufstellen, um die Tiefeninformation lambda fuer die
% 3D Punkte in Cameraframe 1 zu erhalten. Gleiches gilt fuer das GLS M2,
% dass uns die Tiefeninformation fuer die 3D Punkte in CF2 liefert.
lambdas = zeros(N, 2, 4);
posdepth = zeros(1,4);
for i=1:4
    M1 = zeros(3*N,N+1);
    M2 = zeros(3*N,N+1);
    
    for j=1:N
        M1((j-1)*3+1:(j-1)*3+3,j)    = skew(x2_array(:,j)) * R_all(:,:,i) * x1_array(:,j);
        M1((j-1)*3+1:(j-1)*3+3,N+1)  = skew(x2_array(:,j)) * T_all(:,i);
        
        M2((j-1)*3+1:(j-1)*3+3,j)    = skew( R_all(:,:,i) * x1_array(:,j) ) * x2_array(:,j);
        M2((j-1)*3+1:(j-1)*3+3,N+1)  = -skew( R_all(:,:,i) * x1_array(:,j) ) * T_all(:,i);
    end
    
    % Least-Squares-Lösung via Singulärwertzerlegung
    % Da T und damit \gamma nur bis auf Skalierung bestimmbar ist, werden
    % die Tiefen so normiert, das \gamma den Wert 1 erhaelt
    [~,~,V1] = svd(M1,0);
    lambdas1 = V1(:,N+1);
    lambdas1 = lambdas1 ./ lambdas1(end);

    [~,~,V2] = svd(M2,0);
    lambdas2 = V2(:,N+1);
    lambdas2 = lambdas2 ./ lambdas2(end);   

    % Summiere über die Vorzeichen der Tiefen
    posdepth(i) = sum(sign(lambdas1(1:end-1))) + sum(sign(lambdas2(1:end-1)));
    lambdas(:, 1, i)= lambdas1(1:end-1);
    lambdas(:, 2, i)= lambdas2(1:end-1);
end

% Finde die Konfiguration, bei der die meisten Tiefen positiv sind
[~, index] = max(posdepth);

% Wähle das richtige Paar (R,T) und die korrekten Tiefen aus
T       = T_all(:,index);
R       = R_all(:,:,index);

% Zeige euklidische Transformation
fprintf('Rotations-Matrix:\n');
fprintf('%+5.3f %+5.3f %+5.3f\n',R');
fprintf('\nTranslations-Vektor:\n');
fprintf('%+4.2f\n',T);

lambda  = lambdas(:, :, index);

% Vektor mit geschätzten Tiefen aus Ansicht 1
lambda1 = lambda(:,1);

%% Berechnung und Darstellung der 3D-Punkte und der Kameras

% Berechne die Weltkoordinaten der Bildpunkte aus Bild 1 in Kamera 1
P1 = bsxfun(@times, lambda1', x1_array);

% Darstellung der 3D-Rekonstruktion samt Kamera-Perspektiven
figure('name', 'Rekonstruierte 3D-Punkte');

% Zeichne die 3D-Punkte im Kamerasystem 1 ein
hold on
for i = 1:N
    if P1(3,i)>0
        scatter3(P1(1,i), P1(2,i), P1(3,i), '.k');
        text(P1(1,i), P1(2,i), P1(3,i),sprintf('%d', i), 'fontsize', 10, 'color', [0 0 0]);
    end
end
hold off;

% Zeichne beide Kameras ein
camSize = .2;
% Eine Kamera wird durch ein Quadrat im Raum dargestellt
camC1 = [-1 1 1 -1; 1 1 -1 -1]*camSize;
% Homogene Koordinaten der Eckpunkte von Kamera 1
camC1 = [camC1; ones(1,4)];
% Die Koordinaten der zweiten Kamera
camC2 = R\(bsxfun(@plus, camC1, -T));

cam_trace1 = camC1(:,[1:4 1]);
cam_trace2 = camC2(:,[1:4 1]);
hold on;
plot3(cam_trace1(1,:), cam_trace1(2,:), cam_trace1(3,:), 'k');
text(cam_trace1(1,4), cam_trace1(2,4), cam_trace1(3,4), 'C1');
plot3(cam_trace2(1,:), cam_trace2(2,:), cam_trace2(3,:), 'r');
text(cam_trace2(1,4), cam_trace2(2,4), cam_trace2(3,4), 'C2', 'Color', 'r');

xlabel('x'); ylabel('y'); zlabel('z');
grid on; axis equal; hold off;

% Kameraposition und Blickwinkel
campos([43, -22, -87]);
camup([-0.0906, -0.9776, 0.1899]);
camva(3.1874)

end

function [Vhat]=skew(V)
%Umwandlung von V in eine schiefsymmetrische Matrix
Vhat = [0 -V(3) V(2); V(3) 0 -V(1); -V(2) V(1) 0];
end