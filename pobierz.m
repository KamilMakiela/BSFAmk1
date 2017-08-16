function [label, dane, dmu_l, time_l] = pobierz(nvar, n, t)
dane = zeros(n*t,nvar);
i = 0;
while i < nvar
    i = i + 1;
    [mac, dmu_l, time_l, label(i)] = importuj(i);
    if (size(mac,1) ~= n) || (size(mac,2) ~= t)
        i = i - 1;
        h = msgbox('You have chosen a wrong set (inconsistent with n and t values provided earlier). \nReapeat the step. Alternatively you may press ok and ESC to correct n and/or t');
        waitfor(h);
    else
        dane(:,i) = mac(:);
    end
end
%assignin('base', strcat('zmienna',int2str(nvar)), log(dane));
return

function [data, unit_lab, time_lab, var_lab] = importuj(c_nvar)
[FileName,PathName,FilterIndex] = uigetfile({'*.*'},'Open file');
try
    [data,data_label] = xlsread(strcat(PathName, FileName),-1);
catch ex
    disp(ex);
    return;
end
if ~isempty(data_label)
    fprintf('collected data on %d DMUs and %d periods; varable: %s var_num: %d \n', size(data,1), size(data,2), char(data_label(1,1)), c_nvar);
    var_lab = data_label(1,1);
    unit_lab = data_label(2:size(data_label,1),1);
    time_lab = data_label(1,2:size(data_label,2));
else
    fprintf('collected data on %d DMUs and %d perids \nlabels have not been collected!!\n', size(data,1), size(data,2));
    var_lab = cell(1,1);
    unit_lab = cell(2:size(data_label,1),1);
    time_lab = cell(1,2:size(data_label,2));
end
return