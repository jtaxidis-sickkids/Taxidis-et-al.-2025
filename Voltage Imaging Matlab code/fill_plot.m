function fill_plot(x,y,stdy,c)

% Make into rows
if size(x,1) > 1 && size(x,2) == 1                % If not given in row arrays
    x = x';
    y = y';    
    stdy = stdy';
end

% Remove nans
k = isnan(y);
x(k) = [];
y(k) = [];
stdy(k) = [];

% Plot
fill([x, fliplr(x)],[y+stdy, fliplr(y-stdy)],c,'FaceAlpha',0.2,'EdgeColor','none');
hold on;
plot(x,y,'Color',c,'Linewidth',1.5);