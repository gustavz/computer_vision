function [Korrespondenzen] = punkt_korrespondenzen(I1,I2,Mpt1,Mpt2,varargin)
% In dieser Funktion sollen die extrahierten Merkmalspunkte aus einer
% Stereo-Aufnahme mittels NCC verglichen werden um Korrespondenzpunktpaare
% zu ermitteln.

%% Input parser
P = inputParser;

% Liste der optionalen Parameter
% Fensterlänge
P.addOptional('window_length', 25, @isnumeric)
% Minimal geforderte Korrelation
P.addOptional('min_corr', 0.95, @(x)isnumeric(x) && x > 0 && x < 1);
% Plot ein/aus
P.addOptional('do_plot', false, @islogical);

% Lese den Input
P.parse(varargin{:});

% Extrahiere die Variablen aus dem Input-Parser
window_length   = P.Results.window_length;
min_corr        = P.Results.min_corr;
do_plot         = P.Results.do_plot;

Im1             = double(I1);
Im2             = double(I2);

% Es sollen nur die Merkmale berücksichtigt werden, dessen Fenster komplett
% in das Bild hineinpassen. Merkmale, die zu nah am Rand liegen, sollen abgeschnitten werden.
offset_x_begin  = ceil(window_length/2);
offset_x_end    = size(I1,2) - floor(window_length/2);

offset_y_begin  = ceil(window_length/2);
offset_y_end    = size(I1,1) - floor(window_length/2);


% Randbehandlung (lösche alle Merkmale die zu nah am Rand liegen)
% Die Matrix ind1/2 enthaelt in ihren 4 Zeilen an den Stellen eine logische 1, an denen die entsprechende Koordinate zu nah an einem der 4 Bildraender liegt.
ind1            = logical([Mpt1(1,:)<offset_x_begin; Mpt1(1,:)>offset_x_end; Mpt1(2,:)<offset_y_begin; Mpt1(2,:)>offset_y_end]);
ind1            = any(ind1,1);
Mpt1(:,ind1)    = [];   % Loescht die Eintraege an denen die x-, oder die y-Kooridnate zu nah am Rand liegt

ind2            = logical([Mpt2(1,:)<offset_x_begin; Mpt2(1,:)>offset_x_end; Mpt2(2,:)<offset_y_begin; Mpt2(2,:)>offset_y_end]);
ind2            = any(ind2,1);
Mpt2(:,ind2)    = [];   % Loescht die Eintraege an denen die x-, oder die y-Kooridnate zu nah am Rand liegt

% Anzahl der Merkmale in den Bildern
no_pts1         = size(Mpt1,2);
no_pts2         = size(Mpt2,2);

% Zunächst müssen wir die Indizes aller Punkte im Suchfenster bestimmen.
Win_Ix          = -floor(window_length/2):floor(window_length/2);


% Matrix mit allen vektorisierten Merkmalsregionen in Bild 1
Mat_feat_1      = zeros(window_length*window_length,no_pts1);
% Matrix mit allen vektorisierten Merkmalsregionen in Bild 2
Mat_feat_2      = zeros(window_length*window_length,no_pts2);


for index = 1:max(no_pts1,no_pts2)
    if index <= no_pts1
        % Verschiebung der Indizes an korrekte Position
        W_x_coord           = Win_Ix + Mpt1(1,index);
        W_y_coord           = Win_Ix + Mpt1(2,index);
        % Fenster extrahieren
        W_extracted         = Im1(W_y_coord,W_x_coord);
        % Normierung und Abspeichern als Vektor
        Mat_feat_1(:,index) = (W_extracted(:)-mean(W_extracted(:))) / std(W_extracted(:));
    end
    
    if index <= no_pts2
        % Verschiebung der Indizes an korrekte Position
        W_x_coord           = Win_Ix + Mpt2(1,index);
        W_y_coord           = Win_Ix + Mpt2(2,index);
        % Fenster extrahieren
        W_extracted         = Im2(W_y_coord,W_x_coord);
        % Normierung und Abspeichern als Vektor
        Mat_feat_2(:,index) = (W_extracted(:)-mean(W_extracted(:))) / std(W_extracted(:));
    end
end

% Hier wird die letztendliche NCC berechnet durch ein einfaches
% Matrix-Matrix-Produkt. Im Eintrag x,y der NCC_matrix steht die Korrelation des Punktes X
% im zweiten Bild mit dem Punkt Y im ersten Bild
NCC_matrix = 1/(window_length*window_length-1)*Mat_feat_2'*Mat_feat_1;

% Setze alle Korrelationswerte auf 0, die kleiner als der Schwellwert sind
NCC_matrix(NCC_matrix<min_corr) = 0;

% Sortiere die Einträge der Matrix
[sorted_list,sorted_index] = sort(NCC_matrix(:),'descend');

% Eliminiere Elemente aus der Liste, deren Korrelationswert auf null
% gesetzt wurde
sorted_index(sorted_list==0) = [];

% Initialisere die Korrespondenzenmatrix
Korrespondenzen = zeros(4,min(no_pts1,no_pts2));
Korr_count = 1;

% Maximale Anzahl an Iterationen entspricht
max_it      = numel(sorted_index);
size_ncc    = size(NCC_matrix);

for it = 1:max_it    
    % Nehme nächstes Element aus der absteigend nach Größe sortierten Liste
    pt_index = sorted_index(it);
    % Kontrolle ob dieser Wert noch existiert oder bereits auf 0 gesetzt
    % wurde
    if(NCC_matrix(pt_index)==0)
        continue;
    else
        % Extrahiere Reihen- und Spaltenindex
        [Idx_fpt2,Idx_fpt1] = ind2sub(size_ncc,pt_index);
    end
    
    % Setze entsprechende Spalte gleich null, damit ein
    % Merkmalspunkt in Bild 1 nicht mehr als einem anderen Merkmalspunkt in
    % Bild 2 zugewiesen werden kann.
    NCC_matrix(:,Idx_fpt1) = 0;
    
    % Setze entsprechende Zeile gleich null, damit ein
    % Merkmalspunkt in Bild 2 nicht mehr als einem anderen Merkmalspunkt in
    % Bild 1 zugewiesen werden kann.
    % NCC_matrix(Idx_fpt2,:) = 0;
    
    % Speichere das Korrespondenzpunktpaar ab
    Korrespondenzen(:,Korr_count) = [Mpt1(:,Idx_fpt1);Mpt2(:,Idx_fpt2)];
    Korr_count = Korr_count+1;
end
% Lösche überflüssige Elemente der Korrespondenzenmatrix
Korrespondenzen = Korrespondenzen(:,1:Korr_count-1);

%% Zeige die Korrespondenzpunktpaare an
if do_plot
    figure('name', 'Punkt-Korrespondenzen');
    imshow(uint8(I1))
    hold on
    plot(Korrespondenzen(1,:),Korrespondenzen(2,:),'r*')
    imshow(uint8(I2))
    alpha(0.5);
    hold on
    plot(Korrespondenzen(3,:),Korrespondenzen(4,:),'g*')
    for i=1:size(Korrespondenzen,2)
        hold on
        x_1 = [Korrespondenzen(1,i), Korrespondenzen(3,i)];
        x_2 = [Korrespondenzen(2,i), Korrespondenzen(4,i)];
        line(x_1,x_2);
    end
    hold off
end
end


