function SuperTitle(titleText)
% SuperTitle(titleText)
% puts one title at the top of a subplot
% John Grogan, 2019

set(gcf,'NextPlot','add');
axes;
h = title(titleText);
set(gca,'Visible','off');
set(h,'Visible','on');
end