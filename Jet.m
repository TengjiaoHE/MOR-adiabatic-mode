function custom_colormap = Jet()
[all_themes, all_colors] = GetColors();
color1 = all_themes{26};
color2 = all_themes{14};
color3 = all_themes{21};
color = [color1(1,:);color2(1,:);color3(3,:);0.06,0.51,0.04;1,1,0.07;1,0,0;0.7,0,0];
custom_colormap = GenColormap(color, 512);