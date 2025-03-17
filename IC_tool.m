function [IC_data,peaks_infor]=IC_tool(data,options)
% The IC curve analysis tool provides two analysis methods: dVdQ and dQdV.
%-------------------------
% input:
% data: a matrix,size:(N*2),the frist column is capacity,the second colum is
% voltage;
% options is struct with fields:
% method: the analysis method, support 'dVdQ' and 'dQdV', the default is
% 'dQdV';
% point_num: The number of points to interpolate the raw data,size:(1*1), the
% default is 100;
% peaks_plot: Switch to mark the peak point,'off' means not to mark,'on'
% means mark, the default is 'off'.
%-------------------------
%output
% IC_data: the dVdQ and dQdV data,size:(point_num,2);
% peaks_infor: The location of the peak, size(N,2).

arguments
    data (:,2) double ;
    options.method string ='dQdV';
%     options.point_num (1,1) double =100
    options.peaks_plot string='off';
    options.plot_flage string='on'
    options.step (1,1) double=0.005

end
method=options.method;
% point_num=options.point_num;
peaks_plot=options.peaks_plot;
plot_flage= options.plot_flage;
step=options.step ;
switch method
    case 'dQdV'
        x=data(:,2);
        y=data(:,1);
        x_label='Voltage';
        y_label='Capacity';

    case 'dVdQ'
        x=data(:,1);
        y=data(:,2);
        x_label='Capacity';
        y_label='Voltage';
    otherwise
        error('no %s method',method);
end
% x_interp=linspace(min(x),max(x),point_num)';
x_interp= min(x):step:max(x);
x_interp=x_interp(:);
[A,~]=unique(x, 'rows');
if length(A)==length(x)
    y_interp=interp1(x,y,x_interp,"linear","extrap");
else
    [xData, yData] = prepareCurveData( x, y );
    % 设置 fittype 和选项。
    ft = fittype( 'smoothingspline' );
    % 对数据进行模型拟合。
    [fitresult, gof] = fit( xData, yData, ft, 'Normalize', 'on' );
    % 计算插值目标值
    y_interp=feval(fitresult,x_interp);
end

dx=diff(x_interp);
dy=diff(y_interp);

ft = fittype( 'smoothingspline' );
[fitresult_2, ~] = fit( x_interp(2:end),dy./dx, ft, 'Normalize', 'on' );
dydx=feval(fitresult_2,x_interp(2:end));

if data(1,2)<data(end,2)

    switch method
        case 'dVdQ'
            [peaks,peaks_index]=findpeaks(-dydx);
            peaks_infor=[x_interp(1+peaks_index),-peaks];
        case 'dQdV'
            [peaks,peaks_index]=findpeaks(dydx);
            peaks_infor=[x_interp(1+peaks_index),peaks];
    end
else
    switch method
        case 'dVdQ'
            [peaks,peaks_index]=findpeaks(-dydx);
            peaks_infor=[x_interp(1+peaks_index),-peaks];
        case 'dQdV'
            [peaks,peaks_index]=findpeaks(-dydx);
            peaks_infor=[x_interp(1+peaks_index),-peaks];
    end

    %     peaks_infor=[x_interp(1+peaks_index),-peaks];
end
IC_data=[x_interp(2:end),dydx];

% plot
if plot_flage=="on"
    figure
    screenSize = get(0, 'ScreenSize');
    % 计算新的窗口大小
    newWidth = screenSize(3) / 3;
    newHeight = screenSize(4) / 3;
    set(gcf, 'Position', [screenSize(3)/2-newWidth/2 screenSize(4)/2-newHeight/2 newWidth newHeight]);

    subplot(1,2,1)
    plot(x,y,'-o');
    hold on
    plot(x_interp,y_interp,'-*');
    xlabel(x_label);ylabel(y_label);
    grid on;legend({'exp','interp'},Location="northwest")

    subplot(1,2,2)
    plot(x_interp(2:end),dydx,LineWidth=2);
    xlabel(x_label),ylabel(method);grid on;
    hold on
    plot(peaks_infor(:,1),peaks_infor(:,2),'*')
    switch peaks_plot
        case 'on'

            for i = 1:length(peaks_index)
                text(peaks_infor(i,1), peaks_infor(i,2), sprintf('x=%.3f\ny=%.3f', peaks_infor(i,1), peaks_infor(i,2)));
            end
        case 'off'

    end
    % legend({method,'peaks'})
    hold off
else
end
end