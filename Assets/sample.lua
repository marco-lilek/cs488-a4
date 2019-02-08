-- Custom scene

mat1 = gr.material({1, 0.6, 1}, {0.5, 0.7, 0.5}, 25)
mat2 = gr.material({0.6, 1, 0.6}, {0.5, 0.7, 0.5}, 10)
mat3 = gr.material({1.0, 0.6, 0.1}, {0.5, 0.7, 0.5}, 25)
mat4 = gr.material({0.4, 0.4, 0.4}, {1, 1, 1}, 25)
mat5 = gr.material({0.8,0.8,1}, {1, 1, 1}, 25)

scene_root = gr.node('root')
scene_root:translate(0,0,-500)

blayer = gr.node('baselayer')
s1 = gr.sphere('s1')
blayer:add_child(s1)
s1:set_material(mat1)
s1s = 100
s1:scale(s1s, s1s,s1s)
s1:translate(0,-100,0)

s2 = gr.cube('s2')
blayer:add_child(s2)
s2:set_material(mat2)
s2s = 120
s2:scale(s2s,s2s,s2s)
s2:translate(-s2s/2,-s2s/2,-s2s/2)
s2:rotate('z',45)
s2:rotate('y',45)
s2:translate(0,100,0)

helix = gr.node('helix')
scene_root:add_child(helix)

steldodec = gr.mesh( 'dodec', 'smstdodeca.obj' )
steldodec:set_material(mat3)
scene_root:add_child(steldodec)
steldodec:scale(4,4,4)
steldodec:translate(100,-300,-400)

snh = gr.nh_sphere('snh', {500, -300, -400}, 200)
scene_root:add_child(snh)
snh:set_material(mat4)

cnh = gr.nh_box('cnh', {80, -500, -300}, 200)
scene_root:add_child(cnh)
cnh:set_material(mat5)

NUMLAYERS = 10
for i = 1,NUMLAYERS do
  layer = gr.node('layer' .. tostring(i))
  helix:add_child(layer)
  layer:translate(0,0,-100 * i)
  layer:rotate('z',45 * i)

  layer:add_child(blayer)
end
helix:translate(0,0,100*NUMLAYERS/2)
helix:rotate('y',90)
helix:rotate('z',45)
lbboxl = 400.0;
--yellow_light = gr.light({-lbboxl, -lbboxl, lbboxl}, {1, 1, 0}, {1, 0, 0})
red_light = gr.light({-lbboxl, -lbboxl, lbboxl}, {1.0, 0.0, 0.0}, {1, 0, 0})
green_light = gr.light({0, 0, lbboxl}, {0.0, 0.5, 0.0}, {1, 0, 0})
blue_light = gr.light({lbboxl, lbboxl, lbboxl}, {0, 0, 1.0}, {1, 0, 0})

gr.render(scene_root, 'sample.png', 512, 512,
	  {0, 0, 800}, {0, 0, -800}, {0, 1, 0}, 50,
	  {0.3, 0.3, 0.3}, {red_light, blue_light,green_light})
