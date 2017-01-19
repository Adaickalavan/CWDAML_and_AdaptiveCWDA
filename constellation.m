function [S] = constellation(rp)

format = rp.format;
M = rp.M;

%rp.format decides format
switch format
    case 'PSK'% M-PSK, strict Gray coded constellation
        if M == 2 || M == 4 || M == 8 || M == 16
            S(1:M) = exp(1j*2*pi*(0:M-1)/M); %MPSK signal constellation points
            S = S.'; %Transpose to a column vector
        else
            error('Invalid %u-PSK constellation', M)
        end
    case 'QAM'% M-QAM, Rotationally invariant LSB mapping with Gray coding within each sector only
        if M == 16
            %First quadrant
            S(2)=1+1j*3; S(1)=3+1j*3;
            S(4)=1+1j*1; S(3)=3+1j*1;
            %Second quadrant
            S(5)=-3+1j*3; S(7)=-1+1j*3;
            S(6)=-3+1j*1; S(8)=-1+1j*1;
            %Third quadrant
            S(12)=1-1j*1; S(10)=3-1j*1; 
            S(11)=1-1j*3; S(9) =3-1j*3;
            %Fourth quadrant
            S(15)=-3-1j*1; S(16)=-1-1j*1;
            S(13)=-3-1j*3; S(14)=-1-1j*3;
           
            unit=sqrt(1/10);
            S = (S.').*unit; %Transpose to a column vector
        elseif M == 32
            %First quadrant
            S(7)= 1+1j*5; S(3)= 3+1j*5;  
            S(5)= 1+1j*3; S(6)= 3+1j*3; S(8)= 5+1j*3;  
            S(1)= 1+1j*1; S(2)= 3+1j*1; S(4)= 5+1j*1;  
            %Second quadrant
                            S(16)= -3+1j*5; S(12)= -1+1j*5;
            S(11)= -5+1j*3; S(14)= -3+1j*3; S(10)= -1+1j*3;
            S(15)= -5+1j*1; S(13)= -3+1j*1; S(9)= -1+1j*1;
            %Third quadrant
            S(17)= 1-1j*1; S(21)= 3-1j*1; S(23)= 5-1j*1; 
            S(18)= 1-1j*3; S(22)= 3-1j*3; S(19)= 5-1j*3; 
            S(20)= 1-1j*5; S(24)= 3-1j*5; 
            %Fourth quadrant
            S(28)= -5-1j*1; S(26)= -3-1j*1; S(25)= -1-1j*1;
            S(32)= -5-1j*3; S(30)= -3-1j*3; S(29)= -1-1j*3;
                            S(27)= -3-1j*5; S(31)= -1-1j*5;
            
            unit = 1/sqrt(20);
            S = (S.').*unit; %Transpose to a column vector
        elseif M == 64
            %First quadrant
            S(9)= 1+1j*7;  S(10)= 3+1j*7; S(14)= 5+1j*7; S(13)= 7+1j*7;
            S(11)= 1+1j*5; S(12)= 3+1j*5; S(16)= 5+1j*5; S(15)= 7+1j*5;
            S(3)= 1+1j*3;  S(4)= 3+1j*3;  S(8)= 5+1j*3;  S(7)= 7+1j*3;
            S(1)= 1+1j*1;  S(2)= 3+1j*1;  S(6)= 5+1j*1;  S(5)= 7+1j*1;
            %Second quadrant
            S(29)= -7+1j*7; S(31)= -5+1j*7; S(23)= -3+1j*7; S(21)= -1+1j*7;
            S(30)= -7+1j*5; S(32)= -5+1j*5; S(24)= -3+1j*5; S(22)= -1+1j*5;
            S(26)= -7+1j*3; S(28)= -5+1j*3; S(20)= -3+1j*3; S(18)= -1+1j*3;
            S(25)= -7+1j*1; S(27)= -5+1j*1; S(19)= -3+1j*1; S(17)= -1+1j*1;
            %Third quadrant
            S(33)= 1-1j*1; S(35)= 3-1j*1; S(43)= 5-1j*1; S(41)= 7-1j*1;
            S(34)= 1-1j*3; S(36)= 3-1j*3; S(44)= 5-1j*3; S(42)= 7-1j*3;
            S(38)= 1-1j*5; S(40)= 3-1j*5; S(48)= 5-1j*5; S(46)= 7-1j*5;
            S(37)= 1-1j*7; S(39)= 3-1j*7; S(47)= 5-1j*7; S(45)= 7-1j*7;
            %Fourth quadrant
            S(53)= -7-1j*1; S(54)= -5-1j*1; S(50)= -3-1j*1; S(49)= -1-1j*1;
            S(55)= -7-1j*3; S(56)= -5-1j*3; S(52)= -3-1j*3; S(51)= -1-1j*3;
            S(63)= -7-1j*5; S(64)= -5-1j*5; S(60)= -3-1j*5; S(59)= -1-1j*5;
            S(61)= -7-1j*7; S(62)= -5-1j*7; S(58)= -3-1j*7; S(57)= -1-1j*7;
            
            unit = 1/sqrt(42);
            S = (S.').*unit; %Transpose to a column vector
        else
            error('Invalid %u-QAM constellation', M);
        end
    case 'QAMGray' %Strict Gray coded constellation
        if M == 16
            %First quadrant
            S(2)=1+1j*3; S(4)=3+1j*3;
            S(1)=1+1j*1; S(3)=3+1j*1;
            %Second quadrant
            S(8)=-3+1j*3; S(6)=-1+1j*3;
            S(7)=-3+1j*1; S(5)=-1+1j*1;
            %Third quadrant
            S(9)=1-1j*1;  S(11)=3-1j*1; 
            S(10)=1-1j*3; S(12) =3-1j*3;
            %Fourth quadrant
            S(15)=-3-1j*1; S(13)=-1-1j*1;
            S(16)=-3-1j*3; S(14)=-1-1j*3;
           
            unit=sqrt(1/10);
            S = (S.').*unit; %Transpose to a column vector
        elseif M == 64
            %First quadrant
            S(8)= 1+1j*7; S(6)= 3+1j*7; S(14)= 5+1j*7; S(16)= 7+1j*7;
            S(7)= 1+1j*5; S(5)= 3+1j*5; S(13)= 5+1j*5; S(15)= 7+1j*5;
            S(3)= 1+1j*3; S(1)= 3+1j*3; S(9)= 5+1j*3;  S(11)= 7+1j*3;
            S(4)= 1+1j*1; S(2)= 3+1j*1; S(10)= 5+1j*1; S(12)= 7+1j*1;
            %Second quadrant
            S(32)= -7+1j*7; S(30)= -5+1j*7; S(22)= -3+1j*7; S(24)= -1+1j*7;
            S(31)= -7+1j*5; S(29)= -5+1j*5; S(21)= -3+1j*5; S(23)= -1+1j*5;
            S(27)= -7+1j*3; S(25)= -5+1j*3; S(17)= -3+1j*3; S(19)= -1+1j*3;
            S(28)= -7+1j*1; S(26)= -5+1j*1; S(18)= -3+1j*1; S(20)= -1+1j*1;
            %Third quadrant
            S(36)= 1-1j*1; S(34)= 3-1j*1; S(42)= 5-1j*1; S(44)= 7-1j*1;
            S(35)= 1-1j*3; S(33)= 3-1j*3; S(41)= 5-1j*3; S(43)= 7-1j*3;
            S(39)= 1-1j*5; S(37)= 3-1j*5; S(45)= 5-1j*5; S(47)= 7-1j*5;
            S(40)= 1-1j*7; S(38)= 3-1j*7; S(46)= 5-1j*7; S(48)= 7-1j*7;
            %Fourth quadrant
            S(60)= -7-1j*1; S(58)= -5-1j*1; S(50)= -3-1j*1; S(52)= -1-1j*1;
            S(59)= -7-1j*3; S(57)= -5-1j*3; S(49)= -3-1j*3; S(51)= -1-1j*3;
            S(63)= -7-1j*5; S(61)= -5-1j*5; S(53)= -3-1j*5; S(55)= -1-1j*5;
            S(64)= -7-1j*7; S(62)= -5-1j*7; S(54)= -3-1j*7; S(56)= -1-1j*7;

            unit = 1/sqrt(42);
            S = (S.').*unit; %Transpose to a column vector
        else
            error('Invalid %u-QAM Gray constellation', M);
        end
    case 'circular' %Rotationally invariant LSB mapping with Gray code within each sector
        if M == 8 
            %Inner ring
            S(1)= 1+1j;
            S(3)=-1+1j;
            S(5)= 1-1j;
            S(7)=-1-1j;
            %Outer ring
            S(2)=  0+1j*(1+sqrt(3));
            S(4)= -(1+sqrt(3))+1j*0;
            S(6)=  (1+sqrt(3))+1j*0;
            S(8)=  0-1j*(1+sqrt(3));

            unit=1/sqrt(1 + 0.5*(1 + sqrt(3))^2);
            S = (S.').*unit; %Transpose to a column vector
        elseif M == 16
            ring_ratio = 1.59; % Ring ratio
            r1 = sqrt(2/(1 + ring_ratio^2));
            r2 = ring_ratio*r1;
            
            S(1:2) = [r1 r2].*[exp(1j*1*pi/8) exp(1j*0*pi/4)];
            S(3:4) = [r1 r2].*[exp(1j*3*pi/8) exp(1j*1*pi/4)];
            S(7:8) = [r1 r2].*[exp(1j*5*pi/8) exp(1j*2*pi/4)];
            S(5:6) = [r1 r2].*[exp(1j*7*pi/8) exp(1j*3*pi/4)];
            S(13:14) = [r1 r2].*[exp(1j*9*pi/8) exp(1j*4*pi/4)];
            S(15:16) = [r1 r2].*[exp(1j*11*pi/8) exp(1j*5*pi/4)];           
            S(11:12) = [r1 r2].*[exp(1j*13*pi/8) exp(1j*6*pi/4)];
            S(9:10) = [r1 r2].*[exp(1j*15*pi/8) exp(1j*7*pi/4)];
            S = (S.'); %Transpose to a column vector
        else
            error('Invalid %u-cicrcular constellation', M);
        end
    case 'star'
        if M == 16 %Strict Gray coded constellation
            ring_ratio = 2*cosd(67.5)+1; % Ring ratio
            r1 = sqrt(2/(1 + ring_ratio^2));
            r2 = ring_ratio*r1;
            
            S(1:2) = [r1 r2].*exp(1j*0*pi/4);
            S(3:4) = [r1 r2].*exp(1j*1*pi/4);
            S(7:8) = [r1 r2].*exp(1j*2*pi/4);
            S(5:6) = [r1 r2].*exp(1j*3*pi/4);
            S(13:14) = [r1 r2].*exp(1j*4*pi/4);
            S(15:16) = [r1 r2].*exp(1j*5*pi/4);           
            S(11:12) = [r1 r2].*exp(1j*6*pi/4);
            S(9:10) = [r1 r2].*exp(1j*7*pi/4);
            S = (S.'); %Transpose to a column vector
        else
            error('Invalid %u-star constellation', M);
        end
    case 'triangular_NW'
        if M == 8 %No differential encoding available.
            %Inner ring
            S(8)= -2+1j*sqrt(3); S(4)= 0+1j*sqrt(3); S(2)= 2+1j*sqrt(3);         
                             S(3)= -1;         S(1)= 1;
            S(7)= -2-1j*sqrt(3); S(5)= 0-1j*sqrt(3); S(6)= 2-1j*sqrt(3);

            unit = sqrt(2/9);
            S = (S.').*unit; %Transpose to a column vector
        else
            error('Invalid %u-triangular constellation', M);
        end
    case 'triangular_mine'
        if M == 8 %For differential encoding available.
            %Inner ring
            S(3)= -2+1j*sqrt(3); S(4)= 0+1j*sqrt(3); S(2)= 2+1j*sqrt(3);         
                             S(7)= -1;         S(1)= 1;
            S(8)= -2-1j*sqrt(3); S(6)= 0-1j*sqrt(3); S(5)= 2-1j*sqrt(3);

            unit = sqrt(2/9);
            S = (S.').*unit; %Transpose to a column vector
        else
            error('Invalid %u-triangular_mine constellation', M);
        end
    case 'optimum_16'
        %Optimum 16-point constellation from
        %G. J. Foschini, R. D. Gitlin, and S. B. Weinstein, 
        %Optimization of two-dimensional signal constellations in the presence of Gaussian noise
        %IEEE Trans. Commun., vol. 22, no. 1, pp. 28-38, 1974.
        %Table 1, N = 16
        S(1) = 0.007+1i*0.767;
        S(2) = 0.126+1i*0.106;
        S(3) = 0.644+1i*0.545;
        S(4) = 1.279+1i*0.305;
        S(5) = 0.906-1i*0.771;
        S(6) = -1.032-1i*0.103;
        S(7) = -0.504+1i*0.332;
        S(8) = -0.611+1i*1.020;
        S(9) = 0.758-1i*0.119;
        S(10) = -0.911-1i*0.772;
        S(11) = -0.388-1i*0.329;
        S(12) = 0.245-1i*0.552;
        S(13) = -0.272-1i*1.001;
        S(14) = 0.376-1i*1.215;
        S(15) = -1.136+1i*0.571;
        S(16) = 0.512+1i*1.211;
    otherwise
        error('Unavailable %u-%s constellation',M,format);
  
end
    
%------------------------------------------------------------------------------
% %Ensure E[|S|^2] = 1 
% fprintf('signal_Var = %g\n\n', norm(S)^2/M);
% 
% %Get the screensize to specify figure size and location
% scrsz = get(0,'ScreenSize'); 
% %Plot the constellation diagram
% scatterplot(S,1,0,'r*');
% set(gcf,'Outerposition',[2*scrsz(3)/4 scrsz(4)/2 scrsz(3)/4 scrsz(4)/2]);
% title('Signal constellation'),
% xlabel('In-Phase'),ylabel('Qudrature'),grid;  
%------------------------------------------------------------------------------

end

% For 4-PSK (Strict Gray code)
%                             |
%                      (01)->S(2)->sig=1
%                             |
%                             |
% --(11)->S(3)->sig=2---------|--------(00)->S(1)->sig=0--
%                             |
%                             |
%                      (10)->S(4)->sig=3
%                             |

% For 8-PSK (Strict Gray code)
%                             |
%                     3->011->S(3)->sig=2  
%                             |
%       2->010->S(4)->sig=3   |   1->001->S(2)->sig=1
%                             |
% --6->110->S(5)->sig=4-------|--------0->000->S(1)->sig=0--
%                             |
%       7->111->S(6)->sig=5   |   4->100->S(8)->sig=7
%                             |
%                     5->101->S(7)->sig=6  
%                             |

% For 8-triangular (Consists of Gray code penalty. No differential encoding.)
% Bit mapping from:
% Lei Xiao, and Xiaodai Dong
% "The exact transition probability and bit error probability of two-dimensional signalling 
% IEEE Trans Wireless Commun., vol. 4, no. 5, pp. 2600-2609, 2005.
%
%                    |
%                    |
% 7(111)           3(011)          1(001)
%                    |   
%                    |
% ---------2(010)----|----0(000)----------
%                    |
%                    |
% 6(110)           4(100)          5(101)
%                    |
%                    |

% For 8-triangular_mine (Consists of Gray code penalty. Can be tried for differential encoding.)
% Sub-optimum bit to symbol mapping from EE5305 Digital Communication module
%
%                    |
%                    |
% 2(010)           3(011)          1(001)
%                    |   
%                    |
% ---------6(110)----|----0(000)----------
%                    |
%                    |
% 7(111)           5(101)          4(100)
%                    |
%                    |

% For 8-circular (Rotationally invariant LSB mapping with Gray code within each sector)
%                              |
%                           1(001)  
%                              |
%                              |
%                              |
%                              |
%                     2(010)   |   0(000)
%                              |
% --3(011)---------------------|--------------------5(101)--
%                              |
%                     6(110)   |   4(100)
%                              |
%                              |
%                              |
%                              |
%                           7(111)  
%                              |

% For 16-circular (Rotationally invariant LSB mapping with Gray code within each sector)
% Note: This is 16-star constellation with the inner ring rotated by 22.5 degrees
%                                |
%                             7(0111)  
%                                |
%           5(0101)              |              3(0011)
%                                |
%                                |
%                       6(0110)  | 2(0010)
%                                |     
%                    4(0100)     |     0(0000)
% --13(1101)---------------------|--------------------1(0001)--
%                    12(1100)    |     8(1000)
%                                |       
%                       14(1110) | 10(1010)
%                                |
%          15(1111)              |              9(1001)
%                                |
%                                |
%                            11(1011)  
%                                |

% For 16-star (Strict Gray code)
%                                |
%                             7(0111)  
%                                |
%           5(0101)              |              3(0011)
%                                |
%                             6(0110)
%                                |
%                      4(0100)   |    2(0010)
%                                |
% --13(1101)-------12(1100)------|-------0(0000)--------1(0001)--
%                                |
%                     14(1110)   |    8(1000)   
%                                |
%                            10(1010)
%          15(1111)              |              9(1001)
%                                |
%                                |
%                            11(1011)  
%                                |

% For 16-QAMGray (Strict Gray code)
%
%       7(0111)     5(0101)   |    1(0001)    3(0011)
%                             |
%       6(0110)     4(0100)   |    0(0000)    2(0010)
%                             |
%     -----------------------------------------------
%      14(1110)    12(1100)   |    8(1000)   10(1010)   
%                             |
%      15(1111)    13(1101)   |    9(1001)   11(1011)   
%                             |

% For 16-QAM (Rotationally invariant LSB mapping with Gray code within each sector)
%
%      4(0100)     6(0110)    |   1(0001)     0(0000)
%                             |
%      5(0101)     7(0111)    |   3(0011)     2(0010)
%                             |
%     -----------------------------------------------
%      14(1110)    15(1111)   |   11(1011)    9(1001)   
%                             |
%      12(1100)    13(1101)   |   10(1010)    8(1000)   
%                             |

% For 32-QAM (Rotationally invariant LSB mapping with Gray code within each sector)
%                                   |
%             15(01111)  11(01011)  |   6(00110)   2(00010)    
%                                   |
%  10(01010)  13(01101)   9(01001)  |   4(00100)   5(00101)   7(00111)  
%                                   |
%  14(01110)  12(01100)   8(01000)  |   0(00000)   1(00001)   3(00011)   
%                                   |
%  --------------------------------------------------------------------
%                                   |
%  27(11011)  25(11001)  24(11000)  |  16(10000)  20(10100)  22(10110)  
%                                   |
%  31(11111)  29(11101)  28(11100)  |  17(10001)  21(10101)  18(10010)  
%                                   |
%             26(11010)  30(11110)  |  19(10011)  23(10111)    
%                                   |

% For 64-QAMGray (Strict Gray code)
%                                                  |
%  31(011111)  29(011101)  21(010101)  23(010111)  |   7(000111)   5(000101)  13(001101)  15(001111)
%                                                  |
%  30(011110)  28(011100)  20(010100)  22(010110)  |   6(000110)   4(000100)  12(001100)  14(001110)
%                                                  |
%  26(011010)  24(011000)  16(010000)  18(010010)  |   2(000010)   0(000000)   8(001000)  10(001010)
%                                                  |
%  27(011011)  25(011001)  17(010001)  19(010011)  |   3(000011)   1(000001)   9(001001)  11(001011)
%                                                  |
%  -------------------------------------------------------------------------------------------------
%                                                  |
%  59(111011)  57(111001)  49(110001)  51(110011)  |  35(100011)  33(100001)  41(101001)  43(101011)
%                                                  |
%  58(111010)  56(111000)  48(110000)  50(110010)  |  34(100010)  32(100000)  40(101000)  42(101010)
%                                                  |
%  62(111110)  60(111100)  52(110100)  54(110110)  |  38(100110)  36(100100)  44(101100)  46(101110)
%                                                  |
%  63(111111)  61(111101)  53(110101)  55(110111)  |  39(100111)  37(100101)  45(101101)  47(101111)
%                                                  |

% For 64-QAM (Rotationally invariant LSB mapping with Gray code within each sector)
%                                                  |
%  28(011100)  30(011110)  22(010110)  20(010100)  |   8(001000)   9(001001)  13(001101)  12(001100)
%                                                  |
%  29(011101)  31(011111)  23(010111)  21(010101)  |  10(001010)  11(001011)  15(001111)  14(001110)
%                                                  |
%  25(011001)  27(011011)  19(010011)  17(010001)  |   2(000010)   3(000011)   7(000111)   6(000110)
%                                                  |
%  24(011000)  26(011010)  18(010010)  16(010000)  |   0(000000)   1(000001)   5(000101)   4(000100)
%                                                  |
%  -------------------------------------------------------------------------------------------------
%                                                  |
%  52(110100)  53(110101)  49(110001)  48(110000)  |  32(100000)  34(100010)  42(101010)  40(101000)
%                                                  |
%  54(110110)  55(110111)  51(110011)  50(110010)  |  33(100001)  35(100011)  43(101011)  41(101001)
%                                                  |
%  62(111110)  63(111111)  59(111011)  58(111010)  |  37(100101)  39(100111)  47(101111)  45(101101)
%                                                  |
%  60(111100)  61(111101)  57(111001)  56(111000)  |  36(100100)  38(100110)  46(101110)  44(101100)
%                                                  |
