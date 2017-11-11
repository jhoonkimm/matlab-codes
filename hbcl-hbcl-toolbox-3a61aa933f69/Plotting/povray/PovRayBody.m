classdef PovRayBody < RigidBody
  
  properties
    
    povrayIncludeFile = '';
    textureName = 'defaultTexture';
    
  end
  
  methods (Static)
    
    function [] = test()
      %%
      includeFile = 'exampleFiles/Body.inc';
      body = PovRayBody(includeFile);
      body.getRenderingCode();
    end
    
  end
  
  methods
    
    function [this] = PovRayBody(povrayIncludeFileInput)
      %%
      this.povrayIncludeFile = povrayIncludeFileInput;
    end
    
    function [povrayRenderingText] = getRenderingCode(this)
      %%
      [path, name, ext] = fileparts(this.povrayIncludeFile);
      
      povrayRenderingText = [name '(' this.textureName ', '];
      
      R = rodrigues(this.orientation);
      povrayRenderingText = [povrayRenderingText sprintf('%g, ', R(:, 1))];
      povrayRenderingText = [povrayRenderingText sprintf('%g, ', R(:, 2))];
      povrayRenderingText = [povrayRenderingText sprintf('%g, ', R(:, 3))];
      
      povrayRenderingText = [povrayRenderingText sprintf('%g, %g, %g)', this.position)];
    end
    
  end
  
end

