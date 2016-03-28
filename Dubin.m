
line_parameters.vec_r = [0; 1; 0];
line_parameters.vec_q = [1; 1; 0];

n_plot	= 31;
[x_plot, y_plot]= meshgrid(linspace(-10,10,n_plot), linspace(-10,10,n_plot));
line_heading	= acos( line_parameters.vec_q(1) / norm(line_parameters.vec_q) );
commanded_course_heading = zeros(n_plot);

tracking_parameters.course_heading_inf = pi/4;
tracking_parameters.k_path = 0.5;

ell = linspace(-9, 9, 101);
for m1 = 1:101
	line_plot(:, m1) = line_parameters.vec_r + ell(m1)*line_parameters.vec_q;
end

for m1 = 1:n_plot
	for m2 = 1:n_plot
		p_x = x_plot(m1, m2);
		p_y = y_plot(m1, m2);
		
		path_error		= [cos(line_heading) sin(line_heading) 0; ...
			-sin(line_heading) cos(line_heading) 0; 0 0 1]*...
			([p_x; p_y; 0] - line_parameters.vec_r);
		epy				= path_error(2);
		commanded_course_heading(m1, m2) = line_heading - ...
			tracking_parameters.course_heading_inf*(2/pi)*atan(tracking_parameters.k_path * epy);
	end
end


p_x_0	= -10;
p_y_0	= 0;
chi_0	= 90*pi/180;
Vg		= 1;

[t_sim, xi_sim] = ode45(@(t,xi) ode_line_tracking(t, xi, Vg, line_parameters, tracking_parameters), ...
	0:0.1:30, [p_x_0; p_y_0; chi_0; 0]);
for m1 = 1:numel(t_sim)
	cla;
	quiver(x_plot, y_plot, cos(commanded_course_heading), ...
		sin(commanded_course_heading), 0.5, 'MaxHeadSize', 3, 'Color', 'b');
	hold on; axis equal;
	plot(line_plot(1, :), line_plot(2,:), 'LineWidth', 2);
	
	plot(xi_sim(m1, 1), xi_sim(m1, 2), 'k.', 'MarkerSize', 10);
	quiver(xi_sim(m1, 1), xi_sim(m1, 2), cos(xi_sim(m1, 3)), ...
		sin(xi_sim(m1, 3)), 1, 'MaxHeadSize', 3, 'Color', 'k', 'LineWidth', 2);
	drawnow;
end


plot(xi_sim(:, 1), xi_sim(:, 2), 'k.', 'MarkerFaceColor', 'k')
