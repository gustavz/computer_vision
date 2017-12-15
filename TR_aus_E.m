function [T1,R1,T2,R2]=TR_aus_E(E)

% Stelle sicher, dass U und V Rotationsmatrizen sind
[U,S,V]=svd(E);

if det(U) < 0
	U = U*diag([1 1 -1]);
end

if det(V) < 0
	V = V*diag([1 1 -1]);
end

% Der Vektor T liegt im Nullraum von E', ebenso liegt U(:,3) im Nullraum
% von E'. Da die Translation nur bis auf Skalierung geschätzt werden kann,
% können wir diesen Vektor für T verwenden.
T1=U(:,3); 

% T2 zeigt in die entgegengesetzte Richtung
T2=-T1;

RZp=[0 -1 0; 1 0 0;0 0 1];
RZm=[0  1 0;-1 0 0;0 0 1];

R1=U*RZp'*V';
R2=U*RZm'*V';

% Alternativ laesst sich T ueber die Formel aus dem Skript berechnen
T1_hat = U*RZp*S*U';
T1b = [T1_hat(3,2);T1_hat(1,3);T1_hat(2,1);];
T2_hat = U*RZm*S*U';
T2b = [T2_hat(3,2);T2_hat(1,3);T2_hat(2,1)];

end


