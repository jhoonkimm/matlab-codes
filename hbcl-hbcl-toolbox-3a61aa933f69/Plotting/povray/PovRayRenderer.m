classdef PovRayRenderer
  
  properties
    
    cameraRigidBody = RigidBody();
    cameraTargetPoint = [0 0 0];
    displayOutputDuringRender = 1;
    
    povrayRenderingExecutable = '"C:/Program Files/POV-Ray/v3.7/bin/pvengine.exe"';
    
  end
  
  methods (Static)
    
    function [] = test()
      %%
      
      includeFile = 'exampleFiles/Body.inc';
      body = PovRayBody(includeFile);
      body.getRenderingCode();

      renderer = PovRayRenderer();
      
      %       renderer.render(body);
      bodies = {body};
      renderingText = renderer.getPovRayFileText(bodies);
      
      %       fprintf(renderingText);
      
      %       renderer.writePovRayFile('test.pov', bodies);
      
      %       renderer.render(bodies, 'test.png');
      
      imageFilenames = {};
      ts = linspace(0, 2, 30);
      for i = 1 : length(ts)
        t = ts(i);
        body.orientation = body.orientation + [0.1 0.1 -0.05]';
        body.position = [t cos(t), sin(t)]' * 0.01;
        imageFilename = sprintf('exampleFiles/test%g.png', i);
        imageFilenames{i} = imageFilename;
        renderer.render({body}, imageFilename);
      end
      
      frameRate = 15;
      combineImagesIntoVideo(imageFilenames, 'exampleFiles/testVideo2.avi', frameRate);
    end
    
  end

  
  methods
    
    function [this] = PovRayRenderer()
      %%
      this.cameraRigidBody.position = [0.15, 0.15, 0.15];
    end
    
    function [] = render(this, bodies, imageName)
      %%
      [path, name, ext] = fileparts(imageName);

      filename = [path '/' name '.pov'];
      this.writePovRayFile(filename, bodies);
      args = [filename ' +o' imageName];
      
      if (~this.displayOutputDuringRender)
        args = [args ' -D &'];
      end
      
      dos([this.povrayRenderingExecutable ' ' args]);
    end
    
    
    function [] = writePovRayFile(this, filename, bodies)
      %%
      text = this.getPovRayFileText(bodies);
      fid = fopen(filename, 'w');
      fprintf(fid, text);
      fid = fclose(fid);
      %       fid = fclose(fid);
    end
    
    function [text] = getPovRayFileText(this, bodies)
      %%
      text = '';
      text = [text '#include "colors.inc"\n'];
      text = [text '#include "textures.inc"\n'];
      
      for i = 1 : length(bodies)
        body = bodies{i};
        text = [text '#include "' body.povrayIncludeFile '"\n\n'];
      end
      
      location = this.cameraRigidBody.position;
      
      R = rodrigues(this.cameraRigidBody.orientation);
      
      direction = R(:, 1);
      right = R(:, 2);
      up = R(:, 3);
      
      text = [text 'camera {\n'];
      %       text = [text '	location  <0, 0, 0.144600>\n'];
      %       text = [text '	direction <0, 0,  1>\n'];
      %       text = [text '	up        <0, 1,  0>\n'];
      %       text = [text '	right   <4/3, 0,  0>\n'];
      %       text = [text '	look_at   <0, 0, 0>\n'];
      text = [text '	location  <' sprintf('%g, %g, %g', location) '>\n'];
      text = [text '	direction <' sprintf('%g, %g, %g', direction) '>\n'];
      text = [text '	sky <' sprintf('%g, %g, %g', up) '>\n'];
      text = [text '	up        <' sprintf('%g, %g, %g', up) '>\n'];
      %       text = [text '	right   <' sprintf('%g, %g, %g', right) '>\n'];
      text = [text '	look_at   <' sprintf('%g, %g, %g', this.cameraTargetPoint) '>\n'];
      text = [text '}\n\n'];
      
      text = [text 'background { color rgb <1, 1, 1> }\n\n'];
      
      
      %       text = [text 'light_source {<0, 0, 0.144600> color White}\n'];
      %       text = [text 'light_source {<0, 0.144600, 0> color White}\n\n'];
      
      %       text = [text 'light_source {<0, 3, 0> color Blue}\n\n'];
      %       text = [text 'light_source {<3, 0, 0> color Green}\n\n'];
      
      text = [text 'light_source {<0, 3, 0> color White}\n\n'];
      text = [text 'light_source {<3, 0, 0> color White}\n\n'];
      
%       text = [text '#declare defaultTexture = texture {pigment {color red 0.792157 green 0.819608 blue 0.929412} finish {Shiny}}\n\n\n'];
      text = [text '#declare defaultTexture = texture {pigment {color red 0.1 green 0.1 blue 0.1} finish {Shiny}}\n\n\n'];
      
      for i = 1 : length(bodies)
        body = bodies{i};
        text = [text body.getRenderingCode() '\n'];
      end
      
    end
    
  end
  
end

