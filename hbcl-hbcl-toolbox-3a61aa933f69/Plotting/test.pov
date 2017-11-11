#include "colors.inc"
#include "textures.inc"
#include "exampleFiles/Body.inc"

camera {
	location  <0.15, 0.15, 0.15>
	direction <1, 0, 0>
	up        <0, 0, 1>
	right   <0, 1, 0>
	look_at   <0, 0, 0>
}

light_source {<0, 0, 0.144600> color White}
light_source {<0, 0.144600, 0> color White}

#declare defaultTexture = texture {pigment {color red 0.792157 green 0.819608 blue 0.929412} finish {Shiny}}


Body(defaultTexture, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0)
