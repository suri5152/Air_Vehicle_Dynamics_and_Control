function xi_dot = ode_line_tracking(t, xi, Vg, line_parameters, tracking_parameters)

%----- Readability
p_x				= xi(1);
p_y				= xi(2);
course_heading	= xi(3);


%----- Tracking law parameters
course_heading_inf	= tracking_parameters.course_heading_inf;
k_path				= tracking_parameters.k_path;


%----- Closed-loop course heading autopilot characteristics
b1	= 1;
b2	= 1;


%----- Compute lateral error
line_heading	= acos( line_parameters.vec_q(1) / norm(line_parameters.vec_q) );
path_error		= [cos(line_heading) sin(line_heading) 0; ...
	-sin(line_heading) cos(line_heading) 0; 0 0 1]*...
	([p_x; p_y; 0] - line_parameters.vec_r);
epy				= path_error(2);
commanded_course_heading = line_heading - course_heading_inf*(2/pi)*atan(k_path * epy);


xi_dot	= zeros(4, 1);
xi_dot(1, 1)	= Vg*cos(course_heading);
xi_dot(2, 1)	= Vg*sin(course_heading);
xi_dot(3, 1)	= xi(4);
xi_dot(4, 1)	= -b1*xi(4) + b2*(commanded_course_heading - course_heading);