% Full Body local coordinate system creator

clc
clear

% % let user select one c3d file
% 
% [filename, path] = uigetfile('*.c3d', 'Select the .c3d file',...
%     'MultiSelect', 'off');
% 
% % check if uigetfile dialog is not returning zero to prevent error
% if filename ~= 0
% 
%     filepath = cat(2,path,filename);
filepath = 'C:\Users\omar\Aamu M. Sc. TSM\Sciebo\Sina\Data\c3d\Handball\23_Ref_02.c3d';


    tempC3dAcq = btkReadAcquisition(filepath);
    Markers = btkGetMarkers(tempC3dAcq);
    Angles = btkGetAngles(tempC3dAcq);
    Forces = btkGetForces(tempC3dAcq);
    Moments = btkGetMoments(tempC3dAcq);
    Powers = btkGetPowers(tempC3dAcq);
    
    clearvars filepath tempC3dAcq;

%end

clearvars filename path;

%% Definitions

% Calculation of the Sternoclavicular (SC) joint center from the Incisura
% jugularis (IJ) based on the values for a 57 year old right muscular male
% shoulder with an estimated body length of 168 cm provided by 
% Klein Breteler, M. D., Spoor, C. W., & van der Helm, F. C.T. (1999)
% right values were mirrored to the left
USE_KLEIN_BRETELER_SC_JOINT_CENTER_ESTIMATION = true;

% Estimation of the Glenohumeral (GH) joint center from the Incisura jugularis (IJ) 
% based on the values for a 57 year old right muscular male shoulder with an estimated 
% body length of 168 cm provided by Klein Breteler, M. D., Spoor, C. W., & van der Helm, F. C.T. (1999)
% right values were mirrored to the left by changing the sign of the medio lateral axis
USE_KLEIN_BRETELER_GH_JOINT_CENTER_ESTIMATION = false;


%% Lower spine coordinate system

% origin lower spine
O_LS(:,:) = 0.5*(Markers.RPSI(:,:)+Markers.LPSI(:,:));

%cranial axis y
y_axis_tmp_LS (:,:) = Markers.L1(:,:) - O_LS(:,:);

%medio lateral axis z
z_axis_tmp_LS (:,:) = Markers.RPSI(:,:) - O_LS(:,:); 

%anterior posterior axis x
x_axis_LS (:,:) = cross(y_axis_tmp_LS(:,:), z_axis_tmp_LS(:,:));
y_axis_LS (:,:) = cross(z_axis_tmp_LS(:,:), x_axis_LS(:,:));
z_axis_LS (:,:) = cross(x_axis_LS(:,:), y_axis_LS(:,:));

fex_LS = x_axis_LS/norm(x_axis_LS);
fey_LS = y_axis_LS/norm(y_axis_LS);
fez_LS = z_axis_LS/norm(z_axis_LS);
 
%% Upper spine cordinate system

% origin upper spine
O_US(:,:) = (Markers.C7(:,:));

%cranial axis y
y_axis_tmp_US (:,:) = O_US(:,:) - Markers.TH10(:,:);

% %anterior posterior axis x
x_axis_tmp_US(:,:) = (Markers.CLAV(:,:) - O_US(:,:));

%medio lateral axis z
z_axis_US (:,:) = cross(x_axis_tmp_US(:,:),y_axis_tmp_US(:,:)); 

x_axis_US (:,:) = cross(y_axis_tmp_US(:,:),z_axis_US(:,:));

y_axis_US (:,:) = cross(z_axis_US(:,:),x_axis_US(:,:));

fex_US = x_axis_US/norm(x_axis_US);
fey_US = y_axis_US/norm(y_axis_US);
fez_US = z_axis_US/norm(z_axis_US);


%% Right Forearm Coordinate System

O_RFA(:,:) = (Markers.RWRB(:,:));
% + midEMEL(:,:) = 0.5*(Markers.RELBmed(:,:) + Markers.RELB(:,:));

%u1 = y-axis = ISB recommendations: The line connecting US and the midpoint between EL and EM,
%pointing proximally. (Wu, van der Helm et al. 2005)
y_axis_tmp_RFA(:,:) = (0.5*(Markers.RELBmed(:,:) + Markers.RELB(:,:))) - O_RFA(:,:);            
%u2 = x-axis = ISB recommendations: The line perpendicular to the plane through US,RS, 
%and the midpoint between EL and EM, pointing forward.
x_axis_tmp_RFA(:,:) = cross(y_axis_tmp_RFA(:,:),(Markers.RWRA(:,:) - O_RFA(:,:)));            

%u3 = z-axis = ISB recommendations: The common line perpendicular to the
%x_axis and y_axis, pointing to the right
z_axis_RFA(:,:) = cross(x_axis_tmp_RFA(:,:),y_axis_tmp_RFA(:,:));
x_axis_RFA(:,:) = cross(y_axis_tmp_RFA(:,:),z_axis_RFA(:,:));
y_axis_RFA(:,:) = cross(z_axis_RFA(:,:),x_axis_RFA(:,:));

fex_RFA = x_axis_RFA/norm(x_axis_RFA);
fey_RFA = y_axis_RFA/norm(y_axis_RFA);
fez_RFA = z_axis_RFA/norm(z_axis_RFA);

%% Left Forearm Coordinate System

O_LFA(:,:) = (Markers.LWRB(:,:));

%u1 = y-axis = ISB recommendations: The line connecting US and the midpoint between EL and EM,
%pointing proximally. (Wu, van der Helm et al. 2005)
y_axis_tmp_LFA(:,:) = (0.5*(Markers.LELBmed(:,:) + Markers.LELB(:,:))) - O_LFA(:,:);            
%u2 = x-axis = ISB recommendations: The line perpendicular to the plane through US,RS, 
%and the midpoint between EL and EM, pointing forward.
x_axis_tmp_LFA(:,:) = cross(y_axis_tmp_LFA(:,:),(O_LFA(:,:)-Markers.LWRA(:,:)));            

%u3 = z-axis = ISB recommendations: The common line perpendicular to the
%x_axis and y_axis, pointing to the right
z_axis_LFA(:,:) = cross(x_axis_tmp_LFA(:,:),y_axis_tmp_LFA(:,:));
x_axis_LFA(:,:) = cross(y_axis_tmp_LFA(:,:),z_axis_LFA(:,:));
y_axis_LFA(:,:) = cross(z_axis_LFA(:,:),x_axis_LFA(:,:));

fex_LFA = x_axis_LFA/norm(x_axis_LFA);
fey_LFA = y_axis_LFA/norm(y_axis_LFA);
fez_LFA = z_axis_LFA/norm(z_axis_LFA);
%% Thorax Coordinate System

O_T(:,:) = (Markers.CLAV(:,:));
Markers.TH8calc(:,:) = (0.5 * (Markers.TH6(:,:) + Markers.TH10(:,:)));
Markers.TH8STRNMidcalc(:,:) = (0.5 * (Markers.TH8calc(:,:) + Markers.STRN(:,:)));
Markers.C7CLAVMidcalc(:,:) = (0.5 * (Markers.C7(:,:) + Markers.CLAV(:,:)));
y_axis_tmp_T(:,:) = (Markers.C7CLAVMidcalc(:,:) - Markers.TH8STRNMidcalc(:,:));
z_axis_tmp_T(:,:) = cross((O_T(:,:) - Markers.C7CLAVMidcalc(:,:)),...
                          (O_T(:,:) - Markers.TH8STRNMidcalc(:,:)));

x_axis_T(:,:) = cross(y_axis_tmp_T(:,:),z_axis_tmp_T(:,:));
z_axis_T(:,:) = cross(x_axis_T(:,:), y_axis_tmp_T(:,:));
y_axis_T(:,:) = cross(z_axis_T(:,:),x_axis_T(:,:));

fex_T = x_axis_T/norm(x_axis_T);
fey_T = y_axis_T/norm(y_axis_T);
fez_T = z_axis_T/norm(z_axis_T);

%% Right Clavicula Coordinate System

% Calculation of the Sternoclavicular (SC) joint center from the Incisura
% jugularis (IJ) based on the values for a 57 year old right muscular male
% shoulder with an estimated body length of 168 cm provided by 
% Klein Breteler, M. D., Spoor, C. W., & van der Helm, F. C.T. (1999)
if USE_KLEIN_BRETELER_SC_JOINT_CENTER_ESTIMATION == true
    O_RCLAV(:,1) = (Markers.CLAV(:,1)-2.45*10);
    O_RCLAV(:,2) = (Markers.CLAV(:,2)-2.25*10);
    O_RCLAV(:,3) = (Markers.CLAV(:,3)-0.8*10);
else
    O_RCLAV(:,:) = (Markers.CLAV(:,:));
end


%u3 = z-axis = ISB recommendations: The line connecting SC and AC, pointing to AC
z_axis_tmp_RCLAV(:,:) = (Markers.RACR(:,:) - O_RCLAV(:,:));

%u2 = x-axis = ISB recommendations: The line perpendicular to Zc and Yt, pointing
%forward. Note that the Xc-axis is defined with respect to the vertical axis of the thorax (Ytaxis)
%because only two bonylandmarks can be discerned at the clavicle.
x_axis_tmp_RCLAV(:,:) = cross(y_axis_T(:,:),z_axis_tmp_RCLAV(:,:));

%u1 = y-axis = ISB recommendations: The common line perpendicular to the Xc- and Zc-axis, pointing upward.
y_axis_RCLAV(:,:) = cross(z_axis_tmp_RCLAV(:,:),x_axis_tmp_RCLAV(:,:));
x_axis_RCLAV(:,:) = cross(y_axis_RCLAV(:,:),z_axis_tmp_RCLAV(:,:));
z_axis_RCLAV(:,:) = cross(x_axis_RCLAV(:,:), y_axis_RCLAV(:,:));

fex_RCLAV = x_axis_RCLAV/norm(x_axis_RCLAV);
fey_RCLAV = y_axis_RCLAV/norm(y_axis_RCLAV);
fez_RCLAV = z_axis_RCLAV/norm(z_axis_RCLAV);

%% Left Clavicula Coordinate System

% Calculation of the Sternoclavicular (SC) joint center from the Incisura
% jugularis (IJ) based on the values for a 57 year old right muscular male
% shoulder with an estimated body length of 168 cm provided by 
% Klein Breteler, M. D., Spoor, C. W., & van der Helm, F. C.T. (1999)
% right values were mirrored to the left by changing the sign of the 
% medio lateral axis
if USE_KLEIN_BRETELER_SC_JOINT_CENTER_ESTIMATION == true
    O_LCLAV(:,1) = (Markers.CLAV(:,1)- 2.45*10);
    O_LCLAV(:,2) = (Markers.CLAV(:,2)+ 2.25*10);
    O_LCLAV(:,3) = (Markers.CLAV(:,3)- 0.8*10);
else
    O_LCLAV(:,:) = (Markers.CLAV(:,:));
end

%u3 = z-axis = ISB recommendations: The line connecting SC and AC, pointing to AC
z_axis_tmp_LCLAV(:,:) = (Markers.LACR(:,:) - O_LCLAV(:,:));

%u2 = x-axis = ISB recommendations: The line perpendicular to Zc and Yt, pointing
%forward. Note that the Xc-axis is defined with respect to the vertical axis of the thorax (Ytaxis)
%because only two bonylandmarks can be discerned at the clavicle.
x_axis_tmp_LCLAV(:,:) = cross(y_axis_T(:,:),z_axis_tmp_LCLAV(:,:));

%u1 = y-axis = ISB recommendations: The common line perpendicular to the Xc- and Zc-axis, pointing upward.
y_axis_tmp_LCLAV(:,:) = cross(z_axis_tmp_LCLAV(:,:),x_axis_tmp_LCLAV(:,:));
x_axis_LCLAV(:,:) = cross(y_axis_tmp_LCLAV(:,:),z_axis_tmp_LCLAV(:,:));
z_axis_LCLAV(:,:) = cross(x_axis_LCLAV(:,:), y_axis_tmp_LCLAV(:,:));
y_axis_LCLAV(:,:) = cross(z_axis_LCLAV(:,:),x_axis_LCLAV(:,:));

fex_LCLAV = x_axis_LCLAV/norm(x_axis_LCLAV);
fey_LCLAV = y_axis_LCLAV/norm(y_axis_LCLAV);
fez_LCLAV = z_axis_LCLAV/norm(z_axis_LCLAV);

%% Right Shoulder Coordinate System

if USE_KLEIN_BRETELER_GH_JOINT_CENTER_ESTIMATION == true
    O_RSHO(:,1) = (Markers.CLAV(:,1)-8.11*10);
    O_RSHO(:,2) = (Markers.CLAV(:,2)-16.37*10);
    O_RSHO(:,3) = (Markers.CLAV(:,3)-1.79*10);
else
    % mediolateral axis ISB x C3D y
    O_RSHO(:,2) = Markers.RACR(:,2);
    O_RSHO(:,3) = Markers.RSHO(:,3);
    O_RSHO(:,1) = (Markers.RACR(:,1) + Markers.RSHO(:,1))/2;
end

y_axis_tmp_RSHO(:,:) = (O_RSHO(:,:)- y_axis_tmp_RFA(:,:));
z_axis_tmp_RSHO(:,:) = cross(y_axis_tmp_RSHO(:,:),y_axis_RFA);
x_axis_RSHO(:,:) = cross(y_axis_tmp_RSHO(:,:),z_axis_tmp_RSHO(:,:));
y_axis_RSHO(:,:) = cross(z_axis_tmp_RSHO(:,:),x_axis_RSHO(:,:));
z_axis_RSHO(:,:) = cross(x_axis_RSHO(:,:),y_axis_RSHO(:,:));

fex_RSHO = x_axis_RSHO/norm(x_axis_RSHO);
fey_RSHO = y_axis_RSHO/norm(y_axis_RSHO);
fez_RSHO = z_axis_RSHO/norm(z_axis_RSHO);

%% Left Shoulder Coordinate System

if USE_KLEIN_BRETELER_GH_JOINT_CENTER_ESTIMATION == true
    O_LSHO(:,1) = (Markers.CLAV(:,1)-8.11*10);
    O_LSHO(:,2) = (Markers.CLAV(:,2)+16.37*10);
    O_LSHO(:,3) = (Markers.CLAV(:,3)-1.79*10);
else
    % mediolateral axis ISB x C3D y
    O_LSHO(:,2) = Markers.LACR(:,2);
    O_LSHO(:,3) = Markers.LSHO(:,3);
    O_LSHO(:,1) = (Markers.LACR(:,1) + Markers.LSHO(:,1))/2;
end

y_axis_tmp_LSHO(:,:) = (O_LSHO(:,:)- y_axis_tmp_LFA(:,:));
z_axis_tmp_LSHO(:,:) = cross(y_axis_tmp_LSHO(:,:),y_axis_LFA);
x_axis_LSHO(:,:) = cross(y_axis_tmp_LSHO(:,:),z_axis_tmp_LSHO(:,:));
y_axis_LSHO(:,:) = cross(z_axis_tmp_LSHO(:,:),x_axis_LSHO(:,:));
z_axis_LSHO(:,:) = cross(x_axis_LSHO(:,:),y_axis_LSHO(:,:));

fex_LSHO = x_axis_LSHO/norm(x_axis_LSHO);
fey_LSHO = y_axis_LSHO/norm(y_axis_LSHO);
fez_LSHO = z_axis_LSHO/norm(z_axis_LSHO);

%% Right Hand Coordinate System

O_RMANUS(:,:) = Markers.RFIN(:,:);
y_axis_tmp_RMANUS = (0.5*(Markers.RWRA(:,:) + Markers.RWRB(:,:))) - O_RMANUS(:,:);
z_axis_tmp_RMANUS = Markers.RWRA(:,:) - Markers.RWRB(:,:);
x_axis_RMANUS = cross(y_axis_tmp_RMANUS(:,:), z_axis_tmp_RMANUS(:,:));
z_axis_RMANUS = cross(x_axis_RMANUS(:,:), y_axis_tmp_RMANUS(:,:));
y_axis_RMANUS = cross(z_axis_RMANUS(:,:),x_axis_RMANUS(:,:));

fex_RMANUS = x_axis_RMANUS/norm(x_axis_RMANUS);
fey_RMANUS = y_axis_RMANUS/norm(y_axis_RMANUS);
fez_RMANUS = z_axis_RMANUS/norm(z_axis_RMANUS);

%% Left Hand Coordinate System

O_LMANUS(:,:) = Markers.LFIN(:,:);
y_axis_tmp_LMANUS = (0.5*(Markers.LWRA(:,:) + Markers.LWRB(:,:))) - O_LMANUS(:,:);
z_axis_tmp_LMANUS = Markers.LWRA(:,:) - Markers.LWRB(:,:);
x_axis_LMANUS = cross(y_axis_tmp_LMANUS(:,:), z_axis_tmp_LMANUS(:,:));
z_axis_LMANUS = cross(x_axis_LMANUS(:,:), y_axis_tmp_LMANUS(:,:));
y_axis_LMANUS = cross(z_axis_LMANUS(:,:),x_axis_LMANUS(:,:));

fex_LMANUS = x_axis_LMANUS/norm(x_axis_LMANUS);
fey_LMANUS = y_axis_LMANUS/norm(y_axis_LMANUS);
fez_LMANUS = z_axis_LMANUS/norm(z_axis_LMANUS);


%% plots LCS

figure

% upper spine
plot3([O_US(1,1),O_US(1,1)+fex_US(1,1)*1000],...
     [O_US(1,2),O_US(1,2)+fex_US(1,2)*1000],...
     [O_US(1,3),O_US(1,3)+fex_US(1,3)*1000],'r');

hold on

plot3([O_US(1,1),O_US(1,1)+fey_US(1,1)*1000],...
     [O_US(1,2),O_US(1,2)+fey_US(1,2)*1000],...
     [O_US(1,3),O_US(1,3)+fey_US(1,3)*1000],'b');
 
plot3([O_US(1,1),O_US(1,1)+fez_US(1,1)*1000],...
     [O_US(1,2),O_US(1,2)+fez_US(1,2)*1000],...
     [O_US(1,3),O_US(1,3)+fez_US(1,3)*1000],'g');
 
% lower spine

plot3([O_LS(1,1),O_LS(1,1)+fex_LS(1,1)*1000],...
     [O_LS(1,2),O_LS(1,2)+fex_LS(1,2)*1000],...
     [O_LS(1,3),O_LS(1,3)+fex_LS(1,3)*1000],'r');

plot3([O_LS(1,1),O_LS(1,1)+fey_LS(1,1)*1000],...
     [O_LS(1,2),O_LS(1,2)+fey_LS(1,2)*1000],...
     [O_LS(1,3),O_LS(1,3)+fey_LS(1,3)*1000],'b');
 
plot3([O_LS(1,1),O_LS(1,1)+fez_LS(1,1)*1000],...
     [O_LS(1,2),O_LS(1,2)+fez_LS(1,2)*1000],...
     [O_LS(1,3),O_LS(1,3)+fez_LS(1,3)*1000],'g');
 
% thorax

plot3([O_T(1,1),O_T(1,1)+fex_T(1,1)*1000],...
     [O_T(1,2),O_T(1,2)+fex_T(1,2)*1000],...
     [O_T(1,3),O_T(1,3)+fex_T(1,3)*1000],'r');

plot3([O_T(1,1),O_T(1,1)+fey_T(1,1)*1000],...
     [O_T(1,2),O_T(1,2)+fey_T(1,2)*1000],...
     [O_T(1,3),O_T(1,3)+fey_T(1,3)*1000],'b');
 
plot3([O_T(1,1),O_T(1,1)+fez_T(1,1)*1000],...
     [O_T(1,2),O_T(1,2)+fez_T(1,2)*1000],...
     [O_T(1,3),O_T(1,3)+fez_T(1,3)*1000],'g');

% forearm
plot3([O_RFA(1,1),O_RFA(1,1)+fex_RFA(1,1)*1000],...
     [O_RFA(1,2),O_RFA(1,2)+fex_RFA(1,2)*1000],...
     [O_RFA(1,3),O_RFA(1,3)+fex_RFA(1,3)*1000],'r');

plot3([O_RFA(1,1),O_RFA(1,1)+fey_RFA(1,1)*1000],...
     [O_RFA(1,2),O_RFA(1,2)+fey_RFA(1,2)*1000],...
     [O_RFA(1,3),O_RFA(1,3)+fey_RFA(1,3)*1000],'b');
 
plot3([O_RFA(1,1),O_RFA(1,1)+fez_RFA(1,1)*1000],...
     [O_RFA(1,2),O_RFA(1,2)+fez_RFA(1,2)*1000],...
     [O_RFA(1,3),O_RFA(1,3)+fez_RFA(1,3)*1000],'g');
 
 plot3([O_LFA(1,1),O_LFA(1,1)+fex_LFA(1,1)*1000],...
     [O_LFA(1,2),O_LFA(1,2)+fex_LFA(1,2)*1000],...
     [O_LFA(1,3),O_LFA(1,3)+fex_LFA(1,3)*1000],'r');

plot3([O_LFA(1,1),O_LFA(1,1)+fey_LFA(1,1)*1000],...
     [O_LFA(1,2),O_LFA(1,2)+fey_LFA(1,2)*1000],...
     [O_LFA(1,3),O_LFA(1,3)+fey_LFA(1,3)*1000],'b');
 
plot3([O_LFA(1,1),O_LFA(1,1)+fez_LFA(1,1)*1000],...
     [O_LFA(1,2),O_LFA(1,2)+fez_LFA(1,2)*1000],...
     [O_LFA(1,3),O_LFA(1,3)+fez_LFA(1,3)*1000],'g');
 
% right shoulder

plot3([O_RSHO(1,1),O_RSHO(1,1)+fex_RSHO(1,1)*1000],...
     [O_RSHO(1,2),O_RSHO(1,2)+fex_RSHO(1,2)*1000],...
     [O_RSHO(1,3),O_RSHO(1,3)+fex_RSHO(1,3)*1000],'r');

plot3([O_RSHO(1,1),O_RSHO(1,1)+fey_RSHO(1,1)*1000],...
     [O_RSHO(1,2),O_RSHO(1,2)+fey_RSHO(1,2)*1000],...
     [O_RSHO(1,3),O_RSHO(1,3)+fey_RSHO(1,3)*1000],'b');
 
plot3([O_RSHO(1,1),O_RSHO(1,1)+fez_RSHO(1,1)*1000],...
     [O_RSHO(1,2),O_RSHO(1,2)+fez_RSHO(1,2)*1000],...
     [O_RSHO(1,3),O_RSHO(1,3)+fez_RSHO(1,3)*1000],'g');
 
 % left shoulder

plot3([O_LSHO(1,1),O_LSHO(1,1)+fex_LSHO(1,1)*1000],...
     [O_LSHO(1,2),O_LSHO(1,2)+fex_LSHO(1,2)*1000],...
     [O_LSHO(1,3),O_LSHO(1,3)+fex_LSHO(1,3)*1000],'r');

plot3([O_LSHO(1,1),O_LSHO(1,1)+fey_LSHO(1,1)*1000],...
     [O_LSHO(1,2),O_LSHO(1,2)+fey_LSHO(1,2)*1000],...
     [O_LSHO(1,3),O_LSHO(1,3)+fey_LSHO(1,3)*1000],'b');
 
plot3([O_LSHO(1,1),O_LSHO(1,1)+fez_LSHO(1,1)*1000],...
     [O_LSHO(1,2),O_LSHO(1,2)+fez_LSHO(1,2)*1000],...
     [O_LSHO(1,3),O_LSHO(1,3)+fez_LSHO(1,3)*1000],'g');
 
% right hand

plot3([O_RMANUS(1,1),O_RMANUS(1,1)+fex_RMANUS(1,1)*1000],...
     [O_RMANUS(1,2),O_RMANUS(1,2)+fex_RMANUS(1,2)*1000],...
     [O_RMANUS(1,3),O_RMANUS(1,3)+fex_RMANUS(1,3)*1000],'r');

plot3([O_RMANUS(1,1),O_RMANUS(1,1)+fey_RMANUS(1,1)*1000],...
     [O_RMANUS(1,2),O_RMANUS(1,2)+fey_RMANUS(1,2)*1000],...
     [O_RMANUS(1,3),O_RMANUS(1,3)+fey_RMANUS(1,3)*1000],'b');
 
plot3([O_RMANUS(1,1),O_RMANUS(1,1)+fez_RMANUS(1,1)*1000],...
     [O_RMANUS(1,2),O_RMANUS(1,2)+fez_RMANUS(1,2)*1000],...
     [O_RMANUS(1,3),O_RMANUS(1,3)+fez_RMANUS(1,3)*1000],'g');
 
 % Left hand

plot3([O_LMANUS(1,1),O_LMANUS(1,1)+fex_LMANUS(1,1)*1000],...
     [O_LMANUS(1,2),O_LMANUS(1,2)+fex_LMANUS(1,2)*1000],...
     [O_LMANUS(1,3),O_LMANUS(1,3)+fex_LMANUS(1,3)*1000],'r');

plot3([O_LMANUS(1,1),O_LMANUS(1,1)+fey_LMANUS(1,1)*1000],...
     [O_LMANUS(1,2),O_LMANUS(1,2)+fey_LMANUS(1,2)*1000],...
     [O_LMANUS(1,3),O_LMANUS(1,3)+fey_LMANUS(1,3)*1000],'b');
 
plot3([O_LMANUS(1,1),O_LMANUS(1,1)+fez_LMANUS(1,1)*1000],...
     [O_LMANUS(1,2),O_LMANUS(1,2)+fez_LMANUS(1,2)*1000],...
     [O_LMANUS(1,3),O_LMANUS(1,3)+fez_LMANUS(1,3)*1000],'g');
 
 % right clavicle
 
 plot3([O_RCLAV(1,1),O_RCLAV(1,1)+fex_RCLAV(1,1)*1000],...
     [O_RCLAV(1,2),O_RCLAV(1,2)+fex_RCLAV(1,2)*1000],...
     [O_RCLAV(1,3),O_RCLAV(1,3)+fex_RCLAV(1,3)*1000],'r');

plot3([O_RCLAV(1,1),O_RCLAV(1,1)+fey_RCLAV(1,1)*1000],...
     [O_RCLAV(1,2),O_RCLAV(1,2)+fey_RCLAV(1,2)*1000],...
     [O_RCLAV(1,3),O_RCLAV(1,3)+fey_RCLAV(1,3)*1000],'b');
 
plot3([O_RCLAV(1,1),O_RCLAV(1,1)+fez_RCLAV(1,1)*1000],...
     [O_RCLAV(1,2),O_RCLAV(1,2)+fez_RCLAV(1,2)*1000],...
     [O_RCLAV(1,3),O_RCLAV(1,3)+fez_RCLAV(1,3)*1000],'g');
 
 %left clacivle
 
plot3([O_LCLAV(1,1),O_LCLAV(1,1)+fex_LCLAV(1,1)*1000],...
     [O_LCLAV(1,2),O_LCLAV(1,2)+fex_LCLAV(1,2)*1000],...
     [O_LCLAV(1,3),O_LCLAV(1,3)+fex_LCLAV(1,3)*1000],'r');

plot3([O_LCLAV(1,1),O_LCLAV(1,1)+fey_LCLAV(1,1)*1000],...
     [O_LCLAV(1,2),O_LCLAV(1,2)+fey_LCLAV(1,2)*1000],...
     [O_LCLAV(1,3),O_LCLAV(1,3)+fey_LCLAV(1,3)*1000],'b');
 
plot3([O_LCLAV(1,1),O_LCLAV(1,1)+fez_LCLAV(1,1)*1000],...
     [O_LCLAV(1,2),O_LCLAV(1,2)+fez_LCLAV(1,2)*1000],...
     [O_LCLAV(1,3),O_LCLAV(1,3)+fez_LCLAV(1,3)*1000],'g');
 
%% Plot markers
plot3(Markers.RPSI(1,1),Markers.RPSI(1,2),Markers.RPSI(1,3),'o');
plot3(Markers.LPSI(1,1),Markers.LPSI(1,2),Markers.LPSI(1,3),'o');
plot3(Markers.RASI(1,1),Markers.RASI(1,2),Markers.RASI(1,3),'o');
plot3(Markers.LASI(1,1),Markers.LASI(1,2),Markers.LASI(1,3),'o');
plot3(Markers.L1(1,1),Markers.L1(1,2),Markers.L1(1,3),'o');

plot3(Markers.C7(1,1),Markers.C7(1,2),Markers.C7(1,3),'o');
plot3(Markers.CLAV(1,1),Markers.CLAV(1,2),Markers.CLAV(1,3),'o');
plot3(Markers.STRN(1,1),Markers.STRN(1,2),Markers.STRN(1,3),'o');
plot3(Markers.TH6(1,1),Markers.TH6(1,2),Markers.TH6(1,3),'o');
plot3(Markers.TH10(1,1),Markers.TH10(1,2),Markers.TH10(1,3),'o');

plot3(Markers.RANK(1,1),Markers.RANK(1,2),Markers.RANK(1,3),'o');
plot3(Markers.RANKmed(1,1),Markers.RANKmed(1,2),Markers.RANKmed(1,3),'o');

plot3(Markers.LANK(1,1),Markers.LANK(1,2),Markers.LANK(1,3),'o');
plot3(Markers.LANKmed(1,1),Markers.LANKmed(1,2),Markers.LANKmed(1,3),'o');

plot3(Markers.RKNE(1,1),Markers.RKNE(1,2),Markers.RKNE(1,3),'o');
plot3(Markers.RKNEmed(1,1),Markers.RKNEmed(1,2),Markers.RKNEmed(1,3),'o');
plot3(Markers.RTRO(1,1),Markers.RTRO(1,2),Markers.RTRO(1,3),'o');

plot3(Markers.LKNE(1,1),Markers.LKNE(1,2),Markers.LKNE(1,3),'o');
plot3(Markers.LKNEmed(1,1),Markers.LKNEmed(1,2),Markers.LKNEmed(1,3),'o');
plot3(Markers.LTRO(1,1),Markers.LTRO(1,2),Markers.LTRO(1,3),'o');

plot3(Markers.RFH(1,1),Markers.RFH(1,2),Markers.RFH(1,3),'o');
plot3(Markers.LFH(1,1),Markers.LFH(1,2),Markers.LFH(1,3),'o');
plot3(Markers.RBH(1,1),Markers.RBH(1,2),Markers.RBH(1,3),'o');
plot3(Markers.LBH(1,1),Markers.LBH(1,2),Markers.LBH(1,3),'o');

plot3(Markers.RACR(1,1),Markers.RACR(1,2),Markers.RACR(1,3),'o');
plot3(Markers.LACR(1,1),Markers.LACR(1,2),Markers.LACR(1,3),'o');

plot3(Markers.RFRA(1,1),Markers.RFRA(1,2),Markers.RFRA(1,3),'o');
plot3(Markers.LFRA(1,1),Markers.LFRA(1,2),Markers.LFRA(1,3),'o');

plot3(Markers.RUpperArm(1,1),Markers.RUpperArm(1,2),Markers.RUpperArm(1,3),'o');
plot3(Markers.LUpperArm(1,1),Markers.LUpperArm(1,2),Markers.LUpperArm(1,3),'o');

plot3(Markers.RELB(1,1),Markers.RELB(1,2),Markers.RELB(1,3),'o');
plot3(Markers.RELBmed(1,1),Markers.RELBmed(1,2),Markers.RELBmed(1,3),'o');
plot3(Markers.LELB(1,1),Markers.LELB(1,2),Markers.LELB(1,3),'o');
plot3(Markers.LELBmed(1,1),Markers.LELBmed(1,2),Markers.LELBmed(1,3),'o');

plot3(Markers.RWRB(1,1),Markers.RWRB(1,2),Markers.RWRB(1,3),'o');
plot3(Markers.RWRA(1,1),Markers.RWRA(1,2),Markers.RWRA(1,3),'o');
plot3(Markers.LWRB(1,1),Markers.LWRB(1,2),Markers.LWRB(1,3),'o');
plot3(Markers.LWRA(1,1),Markers.LWRA(1,2),Markers.LWRA(1,3),'o');

plot3(Markers.RFRA(1,1),Markers.RFRA(1,2),Markers.RFRA(1,3),'o');
plot3(Markers.LFRA(1,1),Markers.LFRA(1,2),Markers.LFRA(1,3),'o');

plot3(Markers.RSHO(1,1),Markers.RSHO(1,2),Markers.RSHO(1,3),'o');
plot3(Markers.LSHO(1,1),Markers.LSHO(1,2),Markers.LSHO(1,3),'o');

plot3(Markers.RACR(1,1),Markers.RACR(1,2),Markers.RACR(1,3),'o');
plot3(Markers.LACR(1,1),Markers.LACR(1,2),Markers.LACR(1,3),'o');

plot3(Markers.RFIN(1,1),Markers.RFIN(1,2),Markers.RFIN(1,3),'o');
plot3(Markers.LFIN(1,1),Markers.LFIN(1,2),Markers.LFIN(1,3),'o');

%plot3([O_US(1,1),O_LS(1,1)],[O_US(1,2),O_LS(1,2)],[O_US(1,3),O_LS(1,3)],'--k');
%plot3([O_T(1,1),Markers.C7(1,1)],[O_T(1,2),Markers.C7(1,2)],[O_T(1,3),Markers.C7(1,3)],'--k');
%plot3([O_T(1,1),Markers.TH8STRNMidcalc(1,1)],[O_T(1,2),Markers.TH8STRNMidcalc(1,2)],[O_T(1,3),Markers.TH8STRNMidcalc(1,3)],'--k');
%plot3([Markers.C7(1,1),Markers.TH8STRNMidcalc(1,1)],[Markers.C7(1,2),Markers.TH8STRNMidcalc(1,2)],[Markers.C7(1,3),Markers.TH8STRNMidcalc(1,3)],'--k');


%% Plot calculated joint centers and calculated markers

plot3(O_LSHO(1,1),O_LSHO(1,2),O_LSHO(1,3),'+k');
plot3(O_RSHO(1,1),O_RSHO(1,2),O_RSHO(1,3),'+k');

plot3(O_RCLAV(1,1),O_RCLAV(1,2),O_RCLAV(1,3),'+k');
plot3(O_LCLAV(1,1),O_LCLAV(1,2),O_LCLAV(1,3),'+k');

plot3(Markers.TH8STRNMidcalc(1,1),Markers.TH8STRNMidcalc(1,2),Markers.TH8STRNMidcalc(1,3),'ok');
plot3(Markers.C7CLAVMidcalc(1,1),Markers.C7CLAVMidcalc(1,2),Markers.C7CLAVMidcalc(1,3),'ok');
plot3(Markers.TH8calc(1,1),Markers.TH8calc(1,2),Markers.TH8calc(1,3),'ok');


axis equal
grid on
legend('x','y','z');

