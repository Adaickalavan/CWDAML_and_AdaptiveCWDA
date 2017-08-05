%Initilaize differential encoding for 
%8-circular,16-circular,16-star,16-QAMGray,16-QAM,32-QAM,64-QAM,64-QAMGray
%constellations
%%
function [sector_rotation, const_points] = init_DE(S,rp)

    M = rp.M; %Number of signal constellation points
    if (strcmp(rp.format,'triangular_NW') && rp.M == 8) 

        sector_rotation([1 2 4 6]) = 1;
        sector_rotation([3 5 7 8]) = -1;

        const_points([1 2 4 6]) = S([1 2 4 6]); %Quadrant 1
        const_points([3 5 7 8]) = S([3 5 7 8]).*exp(1j*pi); %Quadrant 2
    
    elseif (strcmp(rp.format,'triangular_mine') && rp.M == 8) 

        sector_rotation([1 2 3 4]) = 1;
        sector_rotation([5 6 7 8]) = -1;

        const_points([1 2 3 4]) = S([1 2 3 4]); %Quadrant 1
        const_points([5 6 7 8]) = S([5 6 7 8]).*exp(1j*pi); %Quadrant 2

    elseif (strcmp(rp.format,'circular') && rp.M == 8) || ... 
           (strcmp(rp.format,'QAM') && rp.M == 16) || ... 
           (strcmp(rp.format,'QAMGray') && rp.M == 16) || ...
           (strcmp(rp.format,'QAM') && rp.M == 32) || ...
           (strcmp(rp.format,'QAM') && rp.M == 64) || ... 
           (strcmp(rp.format,'QAMGray') && rp.M == 64)

        L = 4; %Degree of phase ambiguity
        nsp = M/L; %nsp - number of signal points in the sector
        sector_rotation(1:1*nsp) = 1;
        sector_rotation(1*nsp + 1:2*nsp) = 1j;
        sector_rotation(2*nsp + 1:3*nsp) = -1j;
        sector_rotation(3*nsp + 1:4*nsp) = -1;

        const_points(1:1*nsp) = S(1:nsp); %Quadrant 1
        const_points(1*nsp + 1:2*nsp) = S(1*nsp + 1:2*nsp).*exp(-1j*pi/2); %Quadrant 2
        const_points(2*nsp + 1:3*nsp) = S(2*nsp + 1:3*nsp).*exp(1j*pi/2); %Quadrant 3
        const_points(3*nsp + 1:4*nsp) = S(3*nsp + 1:4*nsp).*exp(1j*pi); %Quadrant 4

    elseif (strcmp(rp.format,'star') && rp.M == 16) || ...
           (strcmp(rp.format,'circular') && rp.M == 16)
        
        L = 8;
        nsp = M/L; %nsp - number of signal points in the sector
        sector_rotation(1:1*nsp) = exp(1j*0*pi/4);
        sector_rotation(1*nsp + 1:2*nsp) = exp(1j*1*pi/4);
        sector_rotation(2*nsp + 1:3*nsp) = exp(1j*3*pi/4);
        sector_rotation(3*nsp + 1:4*nsp) = exp(1j*2*pi/4);
        sector_rotation(4*nsp + 1:5*nsp) = exp(1j*7*pi/4);
        sector_rotation(5*nsp + 1:6*nsp) = exp(1j*6*pi/4);
        sector_rotation(6*nsp + 1:7*nsp) = exp(1j*4*pi/4);
        sector_rotation(7*nsp + 1:8*nsp) = exp(1j*5*pi/4);

        const_points(1:1*nsp) = S(1:nsp);
        const_points(1*nsp + 1:2*nsp) = S(1:nsp);
        const_points(2*nsp + 1:3*nsp) = S(1:nsp);
        const_points(3*nsp + 1:4*nsp) = S(1:nsp);
        const_points(4*nsp + 1:5*nsp) = S(1:nsp);
        const_points(5*nsp + 1:6*nsp) = S(1:nsp);
        const_points(6*nsp + 1:7*nsp) = S(1:nsp);
        const_points(7*nsp + 1:8*nsp) = S(1:nsp);
        
    else
        error('Incorrect use of init_DE function');
    end
    
%%
%     %Get the screensize to specify figure size and location
%     scrsz = get(0,'ScreenSize'); 
%     %Plot the constellation diagram
%     quadrant = 4;
%     image = scatterplot(0,1,0,'b.');
%     title('Signal constellation'),
%     xlabel('In-Phase'),ylabel('Qudrature'),grid; 
%     set(gcf,'Outerposition',[1*scrsz(3)/4 scrsz(4)/2 scrsz(3)/4 scrsz(4)/2]);
%     hold on
% %     for i = (quadrant-1)*rp.M/4+(1:rp.M/4)
% %     for i = [1 2 4 8]
%     for i = [3 5 6 7]
%         reply = input('Next ?','S'); %Wait for user input
%         scatterplot(const_points(i),1,0,'r*',image);
%     end
%     hold off

end