function output= three_D_plot_rearrange(input)
n=size(input,1)-1;  %number of rows
m=size(input,2)-1;  %number of columns
x_array=[];         %rows
y_array=[];         %columns
z_array=[];         %data
for i =2:n+1
x_array=[x_array input(1,2:end)];
y_array=[y_array input(i,1)*ones(1,m)];
z_array=[z_array input(i,2:end)];
end
output=[x_array;y_array;z_array];

end