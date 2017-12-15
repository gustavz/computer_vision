function  Merkmale = harris_detektor(Image,varargin) 
% In dieser Funktion soll der Harris-Detektor implementiert werden, der
% Merkmalspunkte aus dem Bild extrahiert


%% Input parser
P = inputParser;

% Liste der optionalen Parameter
% Länge des Bildsegments
P.addOptional('segment_length', 15, @isnumeric)
% Gewichtung für das Harris-Kriterium
P.addOptional('k', 0.05, @isnumeric);
% Schwellwertparameter tau, hängt stark von den Intensitätswerten des
% Bildes ab
P.addOptional('tau', 1e6, @isnumeric);
% Kachelgröße
P.addOptional('tile_size', [200,200], @isnumeric);
% Max. Anzahl Merkmale innerhalb einer Kachel
P.addOptional('N', 5, @isnumeric);
% Minimaler Abstand zwischen zwei Merkmalen
P.addOptional('min_dist', 20, @isnumeric);
% Plot ein/aus
P.addOptional('do_plot', false, @islogical);

% Lese den Input
P.parse(varargin{:});

% Extrahiere die Variablen aus dem Input-Parser
segment_length  = P.Results.segment_length;
k               = P.Results.k;
tau             = P.Results.tau;
tile_size       = P.Results.tile_size;
N               = P.Results.N;
min_dist        = P.Results.min_dist;
do_plot         = P.Results.do_plot;

% Falls bei der Kachelgröße nur die Kantenlänge angegeben wird, verwende
% quadratische Kachel
if numel(tile_size) == 1
    tile_size   = [tile_size,tile_size];
end

%% Vorbereitung zur Feature Detektion
% Approximation des Bildgradienten über das Sobel-Filter
Image           = double(Image);
[Ix,Iy]         = sobel_xy(Image);

% Wir wählen ein Segment mit Gauss-Gewichtung zur mittenbetonten Bestimmung
% der Merkmalsposition. sigma = Filterlänge/5
fmask1          = fspecial('gaussian',[segment_length,1],segment_length/5);

% ZunÃ¤chst werden alle Einträge der Harris-Matrix für jeden Pixel im Bild
% bestimmt. Dies spart viele unnötige Operationen gegenüber einer
% for-Schleife über alle Pixel im Bild!
% G(1,1) ist die Summe aller Ix^2 über das Segment W 
Ixqval          = double(conv2(fmask1,fmask1,Ix.^2, 'same'));
% G(1,1) ist die Summe aller Iy^2 über das Segment W 
Iyqval          = double(conv2(fmask1,fmask1,Iy.^2, 'same'));
% G(1,2) und G(2,1) ist die Summe aller Ix * Iy über das Segment W  
Ixyval          = double(conv2(fmask1,fmask1,Ix.*Iy, 'same'));
 
 
%% Merkmalsextraktion über die Harrismessung
% Harrismessung für alle Pixel des Bildes. (H = det(G) - k*tr(G)^2)
corner          = ((Ixqval.*Iyqval - Ixyval.^2) - k*(Ixqval + Iyqval).^2);

% Bei der vorherigen Faltung wurden die Ränder automatisch mit Nullen
% aufgefüllt, wodurch die Harrismessung im Randbereich des Bildes hohe
% Ausschläge liefert. Diese Werte werden nun mit Null überschrieben.
corner          = corner.*zeroBorder(corner,ceil(segment_length/2));

% Eliminierung von Merkmalen, die kleiner als der Schwellwert tau sind.
corner(corner<=tau)=0;


%% Kontrollmechanismus über den Abstand zweier Merkmale und die maximale Anzahl von Merkmaln in einer Kachel
 
% Die Funktion cake erstellt eine "Kuchenmatrix", die eine kreisförmige
% Anordnung von Nullen beinhaltet und den Rest der Matrix mit Einsen
% auffüllt. Damit können, ausgehend vom stärksten Merkmal, andere Punkte
% unterdrückt werden, die den Mindestabstand hierzu nicht einhalten.
Cake            = cake(min_dist);

% Es muss sichergestellt werden, dass auch die Region um einen Merkmalspunkt am Rand
% elementweise mit der cake-Matrix multipliztiert werden kann. Hierzu muss
% ein Nullrand angefügt werden.
Z               = zeros(size(corner,1)+2*min_dist,size(corner,2)+2*min_dist);
Z((min_dist + 1):(size(corner,1)+min_dist),(min_dist + 1):(size(corner,2)+min_dist)) = corner;
corner          = Z;

% Sortiere alle Merkmale der Stärke nach absteigend
[sorted_list,sorted_index]      = sort(corner(:),'descend');

% Eliminiere alle Einträge, deren Merkmalsstärke auf null gesetzt wurde
sorted_index(sorted_list==0)    = [];

% Anzahl an Merkmalen ungleich null und Größe des Suchfeldes
no_points       = numel(sorted_index);
size_corner     = size(corner);

% Das AKKA ist ein Akkumulatorfeld (Ein Eintrag pro Kachel), welches Aufschluss darüber gibt, wie
% viele Merkmale pro Kachel schon gefunden wurden.
AKKA            = zeros(ceil(size(Image,1)/tile_size(1)),ceil(size(Image,2)/tile_size(2)));

% Alloziere ein Array, in dem die Merkmale gespeichert werden.
Merkmale        = zeros(2,min(numel(AKKA)*N,no_points));
% Feature-Zähler
feature_count   = 1;

for current_point = 1:no_points
    % Nehme nächstes Element aus sortierter Liste
    pt_index    = sorted_index(current_point);
    % Überprüfen, ob dieser Merkmalspunkt noch gültig ist    
    if(corner(pt_index)==0)
        continue;
    else
        % Extrahiere Reihen- und Spalten-Index. Die Matlab-Funktion ind2sub macht das gleiche,
        % benötigt aber mehr Zeit.
        col     = floor(pt_index/size_corner(1));
        row     = pt_index - col*size_corner(1);
        col     = col + 1;
    end

    % Berechnung der Indizes, und damit der ID der zum gefundenen
    % Merkmalspunkt korrespondierenden Kachel Ex und Ey
    Ex          = floor((row-min_dist-1)/(tile_size(1)))+1;
    Ey          = floor((col-min_dist-1)/(tile_size(2)))+1;
    
    % Erhöhe den entsprechenden Eintrag im Akkumulatorarray
    AKKA(Ex,Ey) = AKKA(Ex,Ey)+1;
    
    % Multipliziere Region um den gefundenen Merkmalspunkt elementweise mit der Kuchenmaske
    corner(row-min_dist:row+min_dist,col-min_dist:col+min_dist) = corner(row-min_dist:row+min_dist,col-min_dist:col+min_dist).*Cake;
    
    % Teste, ob die entsprechende Kachel schon genügend Merkmale beinhaltet
    if AKKA(Ex,Ey) == N
        % Falls ja, setzte alle verbleibenden Merkmale innerhalb dieser Kachel auf 0
        corner((((Ex-1)*tile_size(1))+1+min_dist):min(size(corner,1),Ex*tile_size(1)+min_dist),(((Ey-1)*tile_size(2))+1+min_dist):min(size(corner,2),Ey*tile_size(2)+min_dist))=0;   
    end
    
    % Speichere den Merkmalspunkt und berücksichtige dabei den angefügten Nullrand.
    Merkmale(:,feature_count) = [col-min_dist;row-min_dist];
    % Erhöhe den Zähler der Schleife
    feature_count = feature_count+1;
end

% Reduziere die Merkmalsliste auf die gültigen Merkmale
Merkmale = Merkmale(:,1:feature_count-1);

%% Darstellung der gefundenen Merkmale
if do_plot
    figure  
    colormap('gray')
    imagesc(Image)
    hold on;
    plot(Merkmale(1,:), Merkmale(2,:), 'gs');
    plot(Merkmale(1,:), Merkmale(2,:), 'g.');
    axis('off');
end
end

function Mask = zeroBorder(I,W)
    Mask = zeros(size(I));
    Mask((W+1):(size(I,1)-W),(W+1):(size(I,2)-W))=1;
end

function Cake = cake(min_dist)
    [X,Y] = meshgrid(-min_dist:min_dist,[-min_dist:-1,0:min_dist]);
    Cake = sqrt(X.^2+Y.^2)>min_dist;
end
