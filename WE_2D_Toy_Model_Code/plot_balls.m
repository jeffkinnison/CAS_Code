%% Spectral Clustering
load('ball_clustering_15.txt')
ball_clustering = ball_clustering_15;
m = size(ball_clustering, 1);
a = zeros(m, 1);
figure; hold on
for i=1:m
  a(i) = ball_clustering(i, 3);
  plot(ball_clustering(i, 1), ball_clustering(i, 2),'ro', 'MarkerSize', 1.0);
end
s = 0.2; % Marker width in units of X
% Create a scatterplot and return a handle to the 'hggroup' object
h = scatter(ball_clustering(:, 1), ball_clustering(:, 2), 1, a, 'Linewidth', 1.0), colormap(jet), colorbar;
axis([-1 1 -1 1]);
%axis([-1.5 1.5 -0.5 1.25]);
xlabel('x');
ylabel('y');
% Obtain the axes size (in axpos) in points
currentunits = get(gca, 'Units');
set(gca, 'Units', 'Points');
axpos = get(gca, 'Position');
set(gca, 'Units', currentunits);
markerWidth = s/diff(xlim)*axpos(3); % Calculate Marker width in points
set(h, 'SizeData', markerWidth^2);


%% Regular Plot
load('total_weight_on_each_ball_2.txt')
total_weight = total_weight_on_each_ball_2;
m = size(total_weight, 1);
a = zeros(m, 1);
figure; hold on
for i=1:m
  a(i) = total_weight(i, 3);
  plot(total_weight(i, 1),total_weight(i, 2),'ro', 'MarkerSize', 1.0);
end
s = 0.2; % Marker width in units of X
% Create a scatter plot and return a handle to the 'hggroup' object
h = scatter(total_weight(:, 1), total_weight(:, 2), 1, a, 'Linewidth', 1.0), colormap(jet), colorbar;
axis([-1 1 -1 1]);
%axis([-1.5 1.5 -0.5 1.25]);
xlabel('x')
ylabel('y')
% Obtain the axes size (in axpos) in points
currentunits = get(gca, 'Units');
set(gca, 'Units', 'Points');
axpos = get(gca, 'Position');
set(gca, 'Units', currentunits);
markerWidth = s/diff(xlim)*axpos(3); % Calculate Marker width in points
set(h, 'SizeData', markerWidth^2);
