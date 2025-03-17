function [IC_data,peaks_infor]=IC_An(data,options)
% The IC curve analysis tool provides two analysis methods: dVdQ and dQdV.
%-------------------------
% input:
% data: a matrix,size:(N*2),the frist column is capacity,the second colum is
% voltage;
% options is struct with fields:
% method: the analysis method, support 'dVdQ' and 'dQdV', the default is
% 'dQdV';
% point_num: The number of points to interpolate the raw data,size:(1*1), the
% default is 150;
% peaks_plot: Switch to mark the peak point,'off' means not to mark,'on'
% means mark, the default is 'off'.
% minProminence:Minimum protrusion (adjusted according to data
% characteristics),the default is 0;
% windowSize: Smoothing window size (number of data points),the default is
% 10;
% step: differential step size, the  default is 0.01
%-------------------------
%output
% IC_data: the dVdQ and dQdV data,size:(point_num,2);
% peaks_infor: The location of the peak, size(N,2).

arguments
    data (:,2) double ;
    options.method string ='dQdV';
    options.point_num (1,1) double =150
    options.peaks_plot string='off';
    options.plot_flage string='on'
    options.minProminence double =0;
    options.windowSize  (1,1) double=10;
    options.step (1,1) double=0.01
end
method=options.method;
point_num=options.point_num;
peaks_plot=options.peaks_plot;
plot_flage= options.plot_flage;
minProminence= options.minProminence;
windowSize=options.windowSize;
step=options.step;
switch method
    case 'dQdV'
        x=data(:,2);
        y=data(:,1);
        x_label='Voltage';
        y_label='Capacity';

        i=1;
        k=1;
        numV = numel(x);
        while i<numV
            diff=0;
            j=1;
            while diff<step && (i+j)<(numel(y)-1)
                diff=abs(x(i+j)-x(i));
                j=j+1;
            end
            dx(k)=x(i+j)-x(i);
            x_out(k) = x(i);
            dy(k)=y(i+j)-y(i);
            %     i=i+j;
            i=i+1; %try looking at every measurement
            k=k+1;
        end

    case 'dVdQ'
        x=data(:,1);
        y=data(:,2);
        x_label='Capacity';
        y_label='Voltage';

        i=1;
        k=1;
        numV = numel(x);
        while i<numV
            diff=0;
            j=1;
            while diff<step && (i+j)<(numel(y)-1)
                diff=abs(y(i+j)-y(i));
                j=j+1;
            end
            dy(k)=y(i+j)-y(i);
            x_out(k) = x(i);
            dx(k)=x(i+j)-x(i);
            %     i=i+j;
            i=i+1; %try looking at every measurement
            k=k+1;
        end

    otherwise
        error('no %s method',method);
end

% i=1;
% k=1;
% numV = numel(x);
% while i<numV
%     diff=0;
%     j=1;
%     while diff<step && (i+j)<(numel(y)-1)
%         diff=abs(y(i+j)-y(i));
%         j=j+1;
%     end
%     dy(k)=y(i+j)-y(i);
%     x_out(k) = x(i);
%     dx(k)=x(i+j)-x(i);
%     %     i=i+j;
%     i=i+1; %try looking at every measurement
%     k=k+1;
% end

dydx=dy./dx;
% Remove leading 0's and NaN's
NaNindex = find(dydx~=0 & ~isnan(dydx),1);
x_out = x_out(NaNindex:end);
dydx = dydx(NaNindex:end);

%% 过滤-插值
% 原始数据
x_output = x_out;
y_output = dydx;
% point_num = 200;
x_interp = linspace(min(x_output), max(x_output), point_num)';

% 1. 处理重复的 x 点
% 若 x 中有重复点，则对相同 x 的 y 值取平均
[unique_x, ~, ic] = unique(x_output, 'stable');
% 利用 accumarray 对重复的 y 值求平均
unique_y = accumarray(ic, y_output, [], @mean);

%2. 过滤
% 方法一：使用中值滤波平滑数据（medfilt1），对很尖锐的异常进行平滑处理
% windowSize = 10;  % 设置窗口长度，根据数据特点调整
smooth_y = medfilt1(unique_y, windowSize);
% smooth_y = smoothdata(unique_y, 'loess', windowSize);
%卷积滤波
% window = ones(1, windowSize) / windowSize;
% smooth_y = conv(unique_y, window, 'same');

filtered_x = unique_x;
filtered_y = smooth_y;

% 3. 利用平滑样条对处理后的数据进行拟合，计算插值值
[xData, yData] = prepareCurveData(filtered_x(:), filtered_y(:));
ft = fittype('smoothingspline');
[fitresult, gof] = fit(xData, yData, ft, 'Normalize', 'on');

% 计算插值后的目标值
y_interp = feval(fitresult, x_interp);

%% 寻峰
% 设置用于过滤的小峰的参数
% minProminence = 0.2;  % 最小突出度（根据数据特点调整）
if data(1,2)>data(end,2)

    switch method
        case 'dVdQ'
            [peaks,locs]=findpeaks(-y_interp,x_interp, 'MinPeakProminence', minProminence);
            peaks_infor=[locs,-peaks];
        case 'dQdV'
            [peaks,locs]=findpeaks(-y_interp,x_interp, 'MinPeakProminence', minProminence);
            peaks_infor=[locs,-peaks];
    end
else
    switch method
        case 'dVdQ'
            [peaks,locs]=findpeaks(y_interp,x_interp, 'MinPeakProminence', minProminence);
            peaks_infor=[locs,peaks];
        case 'dQdV'
            [peaks,locs]=findpeaks(y_interp,x_interp, 'MinPeakProminence', minProminence);
            peaks_infor=[locs,peaks];
    end

end
IC_data=[x_interp(:),y_interp(:)];


%% 4. 绘图展示
if plot_flage=="on"
    figure
    screenSize = get(0, 'ScreenSize');
    % 计算新的窗口大小
    newWidth = screenSize(3) / 3;
    newHeight = screenSize(4) / 3;
    set(gcf, 'Position', [screenSize(3)/2-newWidth/2 screenSize(4)/2-newHeight/2 newWidth newHeight]);

    % figure
    plot(x_out, dydx, 'b.');           % 显示原始数据（含重复与毛刺）
    hold on; grid on;
    plot(x_interp, y_interp, 'r-', 'LineWidth', 2); % 显示平滑拟合后的曲线

    str_method = char(method);
    xlabel(str_method(end))
    ylabel(method)
    title(sprintf('%s Vs. %s',method,str_method(end)))
    
    switch peaks_plot
        case 'on'
            plot(peaks_infor(:,1), peaks_infor(:,2), 'ko', 'MarkerFaceColor', 'k');
            for i = 1:length(locs)
                text(peaks_infor(i,1), peaks_infor(i,2), sprintf('x=%.3f\ny=%.3f', peaks_infor(i,1), peaks_infor(i,2)));
                legend('原始数据', '平滑拟合曲线','过滤后峰值')
            end
        case 'off'

            legend('原始数据', '平滑拟合曲线');
    end

else

end

end
