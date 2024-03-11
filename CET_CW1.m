clear all;

f = 2.4*10^9; % frequency (given)
c = 3*10^8; % speed of light
beta = 2*pi*f/c; % wave number
n1=1; % Absolute permittivity
n2=3; % Relative permittivity

TX = 0; TY = 4; TZ = 2; % Transmitter [x,y,z] location
RX = 10; RY = 3; RZ = 2; % Receiver [x,y,z] location
movingRX = RX: .01: 22; % Receiver X Point change from point P to Q

% Direct Ray
DD = sqrt((movingRX-TX).^2+(RY-TY)^2); % Direct Distance (Direct Ray from TX to RX)
EDirect = (1./DD).*exp(-j*beta*DD);
EDirect_dB = 20*log10(abs(EDirect)); 

%------------------------------------------------------------------------------------
% First Order Reflection from Top Wall
TX1_Top = 0; TY1_Top = 8;
m1_1 = (TY1_Top-RY)./(TX1_Top-movingRX);
% Find the X reflection point on the wall. Y=mX+C, C=TY1_Top, X=(Y-TY1_Top)/m1_1
RefY1_Top = 6; RefX1_Top = (RefY1_Top-TY1_Top)./m1_1;

D1_1a = sqrt((TX1_Top-movingRX).^2+(TY1_Top-RY)^2); % Straight line distance from TX1_Top to movingPointRX (Sub this in formula)
D1_1b = sqrt((TX-RefX1_Top).^2+(TY-RefY1_Top)^2); % Distance from TX to Reflection Point Top 
D1_1c = sqrt((RefX1_Top-movingRX).^2+(RefY1_Top-RY)^2); % Distance from Reflection Point Top to RX

% Cosine rule to find the Incident angle 
cos_angle1_1 = (D1_1b.^2+D1_1c.^2-DD.^2)./(2*D1_1b.*D1_1c);
cos_angle1_1 = max(-1, min(1, cos_angle1_1));
angleI1_1 = (acosd(cos_angle1_1))./2; % Incident angle
angleT1_1=asind((sind(angleI1_1))./3); % Transmitted angle 
RC1_1 = (n1*cosd(angleI1_1) - n2*cosd(angleT1_1))./(n1*cosd(angleI1_1) + n2*cosd(angleT1_1)); % Reflective coefficient

% Find the electric field loss
E_Ref1_1 = (1./D1_1a).*exp(-j*beta*D1_1a).*RC1_1;
E_Ref1_1_dB = 20*log10(abs(E_Ref1_1));

%------------------------------------------------------------------------------------
% First Order Reflection from Bottom Wall
TX1_Bottom = 0; TY1_Bottom = -4;
m1_2 = (TY1_Bottom-RY)./(TX1_Bottom-movingRX);
% Find the X reflection point on the wall. Y=mX+C, C=TY1_Bottom, X=(Y-TY1_Bottom)/m1_2 
RefY1_Bottom = 0; RefX1_Bottom = (RefY1_Bottom-TY1_Bottom)./m1_2;

D1_2a = sqrt((TX1_Bottom-movingRX).^2+(TY1_Bottom-RY)^2); % Straight line distance from TX1_Bottom to movingPointRX (Sub this in formula)
D1_2b = sqrt((TX-RefX1_Bottom).^2+(TY-RefY1_Bottom)^2); % from TX to Reflection Point Bottom
D1_2c = sqrt((RefX1_Bottom-movingRX).^2+(RefY1_Bottom-RY)^2); % from Reflection Point Bottom to RX

% Cosine rule to find the Incident angle 
cos_angle1_2 = (D1_2b.^2+D1_2c.^2-DD.^2)./(2*D1_2b.*D1_2c);
cos_angle1_2 = max(-1, min(1, cos_angle1_2));
angleI1_2 = (acosd(cos_angle1_2))./2; % Incident angle
angleT1_2=asind((sind(angleI1_2))./3);% Transmitted angle
RC1_2 = (n1*cosd(angleI1_2) - n2*cosd(angleT1_2))./(n1*cosd(angleI1_2) + n2*cosd(angleT1_2));% Reflective coefficient

% Find the electric field loss
E_Ref1_2 = (1./D1_2a).*exp(-j*beta*D1_2a).*RC1_2;
E_Ref1_2_dB = 20*log10(abs(E_Ref1_2));

%------------------------------------------------------------------------------------
% 2nd Order Reflection from Top then Bottom
TX2_Bottom_1 = 0; TY2_Bottom_1 = -8;
m2_1a = (TY2_Bottom_1-RY)./(TX2_Bottom_1-movingRX);
% Find the X reflection point on the wall. Y=mX+C, C=TY2_Bottom_1, X=(Y-TY2_Bottom_1)/m2_1a 
RefY2_Bottom_1 = 0; RefX2_Bottom_1 = (RefY2_Bottom_1-TY2_Bottom_1)./m2_1a;

TX2_Top_1 = 0; TY2_Top_1 = 8;
m2_1b = (TY2_Top_1-RY)./(TX2_Top_1-movingRX);
% Find the X reflection point on the wall. Y=mX+C, C=TY2_Top_1, X=(Y-TY2_Top_1)/m2_1b 
RefY2_Top_1 = 6; RefX2_Top_1 = (RefY2_Top_1-TY2_Top_1)./m2_1b;

D2_1 = sqrt((TX2_Bottom_1-movingRX).^2+(TY2_Bottom_1-RY)^2); % Straight line distance from TX2_Bottom_1 to movingPointRX (Sub this in formula)
D2_1a = sqrt((movingRX-RefX2_Bottom_1).^2+(RY-RefY2_Bottom_1)^2); % from RX to Reflection Point bottom
D2_1b = sqrt((RefX2_Bottom_1-RefX2_Top_1).^2+(RefY2_Bottom_1-RefY2_Top_1)^2); % from Reflection Point bottom to Reflection point top
D2_1c = sqrt((RefX2_Top_1-movingRX).^2+(RefY2_Top_1-RY)^2); % from Reflection Point top to RX
D2_1d = sqrt((TX-RefX2_Top_1).^2+(TY-RefY2_Top_1)^2); % from tx to Reflection Point top
D2_1e = sqrt((TX-RefX2_Bottom_1).^2+(TY-RefY2_Bottom_1)^2); % from tx to Reflection point bottom

% Cosine rule to find the Incident angle
cos_angle2_1 = (D2_1a.^2+D2_1b.^2-D2_1c.^2)./(2*D2_1a.*D2_1b);
cos_angle2_1 = max(-1, min(1, cos_angle2_1));
angleI2_1 = (acosd(cos_angle2_1))./2; % Incident angle
angleT2_1=asind((sind(angleI2_1))./3);% Transmitted angle 
RC2_1 = (n1*cosd(angleI2_1) - n2*cosd(angleT2_1))./(n1*cosd(angleI2_1) + n2*cosd(angleT2_1));% Reflective coefficient

cos_angle2_2 = (D2_1d.^2+D2_1b.^2-D2_1e.^2)./(2*D2_1d.*D2_1b);
cos_angle2_2 = max(-1, min(1, cos_angle2_2));
angleI2_2 = (acosd(cos_angle2_2))./2; % Incident angle
angleT2_2=asind((sind(angleI2_2))./3);% Transmitted angle 
RC2_2 = (n1*cosd(angleI2_2) - n2*cosd(angleT2_2))./(n1*cosd(angleI2_2) + n2*cosd(angleT2_2)); % Reflective coefficient

% Find the electric field loss
E_Ref2_1 = (1./D2_1).*exp(-j*beta*D2_1).*RC2_1.*RC2_2;
E_Ref2_1_dB = 20*log10(abs(E_Ref2_1));

%------------------------------------------------------------------------------------
% 2nd Order Reflection from Bottom then Top
TX2_Top_2 = 0; TY2_Top_2 = 16; 
m2_2a = (TY2_Top_2-RY)./(TX2_Top_2-movingRX);
% Find the X reflection point on the wall. Y=mX+C, C=TY2_Top_2, X=(Y-TY2_Top_2)/m2_2a
RefY2_Top_2 = 6; RefX2_Top_2 = (RefY2_Top_2-TY2_Top_2)./m2_2a;

TX2_Bottom_2 = 0; TY2_Bottom_2 = -4; 
m2_2b = (TY2_Bottom_2-RY)./(TX2_Bottom_2-movingRX);
% Find the X reflection point on the wall. Y=mX+C, C=TY2_Bottom_2, X=(Y-TY2_Bottom_2)/m2_2b
RefY2_Bottom_2 = 0; RefX2_Bottom_2 = (RefY2_Bottom_2-TY2_Bottom_2)./m2_2b;

D2_2 = sqrt((TX2_Top_2-movingRX).^2+(TY2_Top_2-RY)^2); % Straight line distance from TX2_Top to movingPointRX (Sub this in formula)
D2_2a = sqrt((movingRX-RefX2_Top_2).^2+(RY-RefY2_Top_2)^2); % from RX to Reflection Point top 
D2_2b = sqrt((RefX2_Bottom_2-RefX2_Top_2).^2+(RefY2_Bottom_2-RefY2_Top_2)^2); % from Reflection Point bottom to Reflection point top
D2_2c = sqrt((RefX2_Bottom_2-movingRX).^2+(RefY2_Bottom_2-RY)^2); % from Reflection Point bottom to RX
D2_2d = sqrt((TX-RefX2_Top_2).^2+(TY-RefY2_Top_2)^2); % from tx to Reflection Point top
D2_2e = sqrt((TX-RefX2_Bottom_2).^2+(TY-RefY2_Bottom_2)^2); % from tx to Reflection point bottom

% Cosine rule to find the Incident angle
cos_angle2_a = (D2_2a.^2+D2_2b.^2-D2_2c.^2)./(2*D2_2a.*D2_2b);
cos_angle2_a = max(-1, min(1, cos_angle2_a));
angleI2_a = (acosd(cos_angle2_a))./2; % Incident angle
angleT2_a=asind((sind(angleI2_a))./3);% Transmitted angle 
RC2_2a = (n1*cosd(angleI2_a) - n2*cosd(angleT2_a))./(n1*cosd(angleI2_a) + n2*cosd(angleT2_a));% Reflective coefficient

cos_angle2_b = (D2_2b.^2+D2_2d.^2-D2_2e.^2)./(2*D2_2b.*D2_2d);
cos_angle2_b = max(-1, min(1, cos_angle2_b));
angleI2_b = (acosd(cos_angle2_b))./2; % Incident angle
angleT2_b=asind((sind(angleI2_b))./3); % Transmitted angle 
RC2_2b = (n1*cosd(angleI2_b) - n2*cosd(angleT2_b))./(n1*cosd(angleI2_b) + n2*cosd(angleT2_b));% Reflective coefficient

% Find the electric field loss
E_Ref2_2 = (1./D2_2).*exp(-j*beta*D2_2).*RC2_2a.*RC2_2b;
E_Ref2_2_dB = 20*log10(abs(E_Ref2_2));

%------------------------------------------------------------------------------------
% 3rd Order Reflection from top then bottom then top
TX3_Top_1a = 0; TY3_Top_1a = 20;
m3_1a = (TY3_Top_1a-RY)./(TX3_Top_1a-movingRX);
% Find the X reflection point on the wall. Y=mX+C, C=TY1_Top_1, X=(Y-TY3_Top_1a)/m3_1a
RefY3_Top_1a = 6; RefX3_Top_1a = (RefY3_Top_1a-TY3_Top_1a)./m3_1a;

TX3_Bottom_1 = 0; TY3_Bottom_1 = -8;
m3_1b = (TY3_Bottom_1-RefY3_Top_1a)./(TX3_Bottom_1-RefX3_Top_1a);
% Find the X reflection point on the wall. Y=mX+C, C=TY3_Bottom_1, X=(Y-TY3_Bottom_1)/m3_1b
RefY3_Bottom_1 = 0; RefX3_Bottom_1 = (RefY3_Bottom_1-TY3_Bottom_1)./m3_1b;

TX3_Top_1b = 0; TY3_Top_1b = 8;
m3_1c = (TY3_Top_1b-RefY3_Bottom_1)./(TX3_Top_1b-RefX3_Bottom_1);
% Find the X reflection point on the wall. Y=mX+C, C=TY1_Top_1b, X=(Y-TY3_Top_1b)/m3_1c
RefY3_Top_1b = 6; RefX3_Top_1b = (RefY3_Top_1b-TY3_Top_1b)./m3_1c;

D3_1 = sqrt((TX3_Top_1a-movingRX).^2+(TY3_Top_1a-RY)^2); % Straight line distance from TX3_Top_1a to movingPointRX (Sub this in formula)
D3_1a = sqrt((movingRX-RefX3_Top_1a).^2+(RY-RefY3_Top_1a)^2); % from RX to Reflection Point Top 1a
D3_1b = sqrt((RefX3_Top_1a-RefX3_Bottom_1).^2+(RefY3_Top_1a-RefY3_Bottom_1 )^2); % from Reflection Point bottom to Reflection point top 1a
D3_1c = sqrt((RefX3_Bottom_1-movingRX).^2+(RefY3_Bottom_1-RY)^2); % from Reflection Point Bottom to RX
D3_1d = sqrt((RefX3_Top_1b -RefX3_Top_1a).^2+(RefY3_Top_1b-RefY3_Top_1a)^2); % from Reflection point Top 1a to Reflection Point Top 1b
D3_1e = sqrt((RefX3_Top_1b-RefX3_Bottom_1).^2+(RefY3_Top_1b-RefY3_Bottom_1 )^2); % from Reflection Point bottom to Reflection point top 1b
D3_1f = sqrt((RefX3_Bottom_1-TX).^2+(RefY3_Bottom_1-TY)^2); % from Reflection Point Bottom to TX)
D3_1g = sqrt((TX-RefX3_Top_1b).^2+(TY-RefY3_Top_1b)^2); % from TX to Reflection Point Top 1b

% Cosine rule to find the Incident angle
cos_angle3_1a = (D3_1a.^2+D3_1b.^2-D3_1c.^2)./(2*D3_1a.*D3_1b);
cos_angle3_1a = max(-1, min(1, cos_angle3_1a));
angleI3_1a = (acosd(cos_angle3_1a))./2; % Incident angle
angleT3_1a = asind((sind(angleI3_1a))./3);% Transmitted angle 
RC3_1a = (n1*cosd(angleI3_1a) - n2*cosd(angleT3_1a))./(n1*cosd(angleI3_1a) + n2*cosd(angleT3_1a));% Reflective coefficient

cos_angle3_1b = (D3_1b.^2+D3_1e.^2-D3_1d.^2)./(2*D3_1b.*D3_1e);
cos_angle3_1b = max(-1, min(1, cos_angle3_1b));
angleI3_1b = (acosd(cos_angle3_1b))./2; % Incident angle
angleT3_1b = asind((sind(angleI3_1b))./3);% Transmitted angle 
RC3_1b = (n1*cosd(angleI3_1b) - n2*cosd(angleT3_1b))./(n1*cosd(angleI3_1b) + n2*cosd(angleT3_1b));% Reflective coefficient

cos_angle3_1c = (D3_1e.^2+D3_1g.^2-D3_1f.^2)./(2*D3_1e.*D3_1g);
cos_angle3_1c = max(-1, min(1, cos_angle3_1c));
angleI3_1c = (acosd(cos_angle3_1c))./2; % Incident angle
angleT3_1c = asind((sind(angleI3_1c))./3);% Transmitted angle 
RC3_1c = (n1*cosd(angleI3_1c) - n2*cosd(angleT3_1c))./(n1*cosd(angleI3_1c) + n2*cosd(angleT3_1c));% Reflective coefficient

% Find the electric field loss
E_Ref3_1 = (1./D3_1).*exp(-j*beta*D3_1).*RC3_1a.*RC3_1b.*RC3_1c;
E_Ref3_1_dB = 20*log10(abs(E_Ref3_1));

%------------------------------------------------------------------------------------
% 3rd Order Reflection from bottom then top then bottom

TX3_Bottom_2a = 0; TY3_Bottom_2a = -16;
m3_2a = (TY3_Bottom_2a-RY)./(TX3_Bottom_2a-movingRX);
% Find the X reflection point on the wall. Y=mX+C, C=TY3_Bottom_2a, X=(Y-TY3_Bottom_2a)/m3_2a
RefY3_Bottom_2a = 0; RefX3_Bottom_2a = (RefY3_Bottom_2a-TY3_Bottom_2a)./m3_2a;

TX3_Top_2 = 0; TY3_Top_2 = 16;
m3_2b = (TY3_Top_2-RefY3_Bottom_2a)./(TX3_Top_2-RefX3_Bottom_2a);
% Find the X reflection point on the wall. Y=mX+C, C=TY3_Top_2a, X=(Y-TY3_Top_2a)/m3_2b
RefY3_Top_2 = 0; RefX3_Top_2 = (RefY3_Top_2-TY3_Top_2)./m3_2b;

TX3_Bottom_2b = 0; TY3_Bottom_2b = -4;
m3_2c = (TY3_Bottom_2b-RefY3_Top_2)./(TX3_Bottom_2b-RefX3_Top_2);
% Find the X reflection point on the wall. Y=mX+C, C=TY3_Bottom_2b, X=(Y-TY3_Bottom_2b)/m3_2c
RefY3_Bottom_2b = 6; RefX3_Bottom_2b = (RefY3_Bottom_2b-TY3_Bottom_2b)./m3_2c;

D3_2 = sqrt((TX3_Bottom_2a-movingRX).^2+(TY3_Bottom_2a-RY)^2); % Straight line distance from TX3_Bottom_2a to movingPointRX (Sub this in formula)
D3_2a = sqrt((movingRX-RefX3_Bottom_2a).^2+(RY-RefY3_Bottom_2a)^2); % from RX to Reflection Point bottom 2a
D3_2b = sqrt((RefX3_Top_2-RefX3_Bottom_2a).^2+(RefY3_Top_2-RefY3_Bottom_2a )^2); % from Reflection Point bottom 2a to Reflection point top
D3_2c = sqrt((RefX3_Top_2-movingRX).^2+(RefY3_Top_2-RY)^2); % from Reflection Point top to RX
D3_2d = sqrt((RefX3_Bottom_2b -RefX3_Bottom_2a).^2+(RefY3_Bottom_2b-RefY3_Bottom_2a)^2); % from Reflection Point bottom 2a to Reflection Point bottom 2b
D3_2e = sqrt((RefX3_Top_2-RefX3_Bottom_2b).^2+(RefY3_Top_2-RefY3_Bottom_2b )^2); % from Reflection Point bottom 2b to Reflection point top
D3_2f = sqrt((RefX3_Top_2-TX).^2+(RefY3_Top_2-TY)^2); % from Reflection Point top to TX
D3_2g = sqrt((TX-RefX3_Bottom_2b).^2+(TY-RefY3_Bottom_2b)^2); % from TX to Reflection Point bottom 2b

% Cosine rule to find the Incident angle
cos_angle3_2a = (D3_2a.^2+D3_2b.^2-D3_2c.^2)./(2*D3_2a.*D3_2b);
cos_angle3_2a = max(-1, min(1, cos_angle3_2a));
angleI3_2a = (acosd(cos_angle3_2a))./2; % Incident angle
angleT3_2a = asind((sind(angleI3_2a))./3);% Transmitted angle 
RC3_2a = (n1*cosd(angleI3_2a) - n2*cosd(angleT3_2a))./(n1*cosd(angleI3_2a) + n2*cosd(angleT3_2a));% Reflective coefficient

cos_angle3_2b = (D3_2b.^2+D3_2e.^2-D3_2d.^2)./(2*D3_2b.*D3_2e);
cos_angle3_2b = max(-1, min(1, cos_angle3_2b));
angleI3_2b = (acosd(cos_angle3_2b))./2; % Incident angle
angleT3_2b = asind((sind(angleI3_2b))./3);% Transmitted angle 
RC3_2b = (n1*cosd(angleI3_2b) - n2*cosd(angleT3_2b))./(n1*cosd(angleI3_2b) + n2*cosd(angleT3_2b));% Reflective coefficient

cos_angle3_2c = (D3_2e.^2+D3_2g.^2-D3_2f.^2)./(2*D3_2e.*D3_2g);
cos_angle3_2c = max(-1, min(1, cos_angle3_2c));
angleI3_2c = (acosd(cos_angle3_2c))./2; % Incident angle
angleT3_2c = asind((sind(angleI3_2c))./3);% Transmitted angle 
RC3_2c = (n1*cosd(angleI3_2c) - n2*cosd(angleT3_2c))./(n1*cosd(angleI3_2c) + n2*cosd(angleT3_2c));% Reflective coefficient

% Find the electric field loss
E_Ref3_2 = (1./D3_2).*exp(-j*beta*D3_2).*RC3_2a.*RC3_2b.*RC3_2c;
E_Ref3_2_dB = 20*log10(abs(E_Ref3_2));

%------------------------------------------------------------------------------------
% Calculate total electric field 
% Direct + 1st Order
Total1st=E_Ref1_2+E_Ref1_1+EDirect;
Total1st_dB= 20*log10(abs(Total1st));
% Direct + 1st Order + 2nd Order
Total2nd=EDirect+E_Ref1_1+E_Ref1_2+E_Ref2_1+E_Ref2_2;
Total2nd_dB= 20*log10(abs(Total2nd));
% Direct + 1st Order + 2nd Order + 3rd Order
Total3rd=EDirect+E_Ref1_1+E_Ref1_2+E_Ref2_1+E_Ref2_2+E_Ref3_1+E_Ref3_2;
Total3rd_dB= 20*log10(abs(Total3rd));

% Plot the graph
plot(movingRX, EDirect_dB, 'LineWidth', 1.5, 'DisplayName', 'Direct Ray');
hold on;
plot(movingRX, E_Ref1_1_dB, 'DisplayName', '1st Order Reflection from Top Wall');
plot(movingRX, E_Ref1_2_dB, 'DisplayName', '1st Order Reflection from Bottom Wall');
plot(movingRX, E_Ref2_1_dB, 'DisplayName', '2nd Order Reflection from Top then Bottom');
plot(movingRX, E_Ref2_2_dB, 'DisplayName', '2nd Order Reflection from Bottom then Top');
plot(movingRX, E_Ref3_1_dB, 'DisplayName', '3rd Order Reflection from Top then Bottom then Top');
plot(movingRX, E_Ref3_2_dB, 'DisplayName', '3rd Order Reflection from Bottom then Top then Bottom');
plot(movingRX, Total1st_dB, 'DisplayName', 'Direct + 1st Order');
plot(movingRX, Total2nd_dB, 'DisplayName', 'Direct + 1st Order + 2nd Order');
plot(movingRX, Total3rd_dB, 'DisplayName', 'Direct + 1st Order + 2nd Order + 3rd Order');

legend('show');
xlabel('Distance (from point P to point Q)');
ylabel('Electric field strength (dB)');
title('Graph of Electric field strength (dB) against Distance (from point P to point Q)' );
grid on;




