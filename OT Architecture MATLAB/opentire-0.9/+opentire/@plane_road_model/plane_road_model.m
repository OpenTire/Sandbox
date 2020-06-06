%PLANE_ROAD_MODEL Simple plane road surface model
classdef (Sealed = true) plane_road_model < opentire.road_model

properties
origin = [0;0;0];
normal = [0;0;1];
grip = 1;
end

methods
  
% compute_contact - compute road-tire contact point, normal, curvature, transforms, and grip
% W(:,3) is road normal
% R(:,2) is wheel_axis
function [tire_contact, road_curvature, R, W, road_grip] = compute_contact(road, wheel_axis, wheel_center, time, distance)
  
  road_normal = road.normal;

  R = zeros(3,3);
  R(:,1) = vec3_normalize(vec3_cross(wheel_axis, road_normal));
  R(:,2) = wheel_axis;
  R(:,3) = vec3_cross(R(:,1), R(:,2));
  
  % W = transform from Tydex W tire frame to vehicle frame
  % W(:,3) == road surface normal, pointing upwards
  % 
  W = zeros(3,3);
  % Z axis: road surface normal, pointing upwards.
  W(:,3) = road_normal;
  % X axis: intersection of wheel plane with ground plane
  W(:,1) = vec3_normalize(vec3_cross(wheel_axis, road_normal));
  % Y axis: projection of wheel axis onto ground plane
  W(:,2) = vec3_cross(road_normal, W(:,1));

  loaded_radius = vec3_dot(wheel_center - road.origin, road.normal) / vec3_dot(R(:,3), road.normal);
  radial = -loaded_radius * R(:,3); % vector from center to contact, in car frame
  tire_contact = wheel_center + radial;
  
  road_grip = road.grip;
  road_curvature = [0;0];
end

end % methods

end % classdef

function v = vec3_normalize(v)
  n = sqrt(v(1)^2 + v(2)^2 + v(3)^2);
  % norm assumed to be > 0
  %assert(n > 0);
  v = v / n;
end

function d = vec3_dot(a, b)
  d = a(1)*b(1) + a(2)*b(2) + a(3)*b(3);
end

function c = vec3_cross(a, b)
  c = [a(2)*b(3) - a(3)*b(2);
       a(3)*b(1) - a(1)*b(3);
       a(1)*b(2) - a(2)*b(1)];
end

