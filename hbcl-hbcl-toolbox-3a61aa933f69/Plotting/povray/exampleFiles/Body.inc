

#ifndef( testPart_Inc_Temp)
#declare testPart_Inc_Temp = version;
#version 3.6;

#macro Body(textureToUse, Val00, Val01, Val02,
        Val10, Val11, Val12,
        Val20, Val21, Val22,
        Val30, Val31, Val32)
object {
  // _X1_X2
  union {
       triangle {
      <-0.024100, 0.005000, 0.021762>,
      <-0.024100, 0.005000, -0.021762>,
      < -0.024100, -0.005000, 0.021762>
    }
    triangle {
      <-0.024100, -0.005000, 0.021762>,
      <-0.024100, 0.005000, -0.021762>,
      < -0.024100, -0.005000, -0.021762>
    }
    triangle {
      <0.024100, 0.005000, 0.021762>,
      <-0.024100, 0.005000, 0.021762>,
      < 0.024100, -0.005000, 0.021762>
    }
    triangle {
      <0.024100, -0.005000, 0.021762>,
      <-0.024100, 0.005000, 0.021762>,
      < -0.024100, -0.005000, 0.021762>
    }
    triangle {
      <0.024100, 0.005000, -0.021762>,
      <0.024100, 0.005000, 0.021762>,
      < 0.024100, -0.005000, -0.021762>
    }
    triangle {
      <0.024100, -0.005000, -0.021762>,
      <0.024100, 0.005000, 0.021762>,
      < 0.024100, -0.005000, 0.021762>
    }
    triangle {
      <-0.024100, 0.005000, -0.021762>,
      <0.024100, 0.005000, -0.021762>,
      < -0.024100, -0.005000, -0.021762>
    }
    triangle {
      <-0.024100, -0.005000, -0.021762>,
      <0.024100, 0.005000, -0.021762>,
      < 0.024100, -0.005000, -0.021762>
    }
    triangle {
      <0.024100, 0.005000, 0.021762>,
      <0.024100, 0.005000, -0.021762>,
      < -0.024100, 0.005000, 0.021762>
    }
    triangle {
      <-0.024100, 0.005000, 0.021762>,
      <0.024100, 0.005000, -0.021762>,
      < -0.024100, 0.005000, -0.021762>
    }
    triangle {
      <0.024100, -0.005000, -0.021762>,
      <0.024100, -0.005000, 0.021762>,
      < -0.024100, -0.005000, -0.021762>
    }
    triangle {
      <-0.024100, -0.005000, -0.021762>,
      <0.024100, -0.005000, 0.021762>,
      < -0.024100, -0.005000, 0.021762>
    }
    texture{textureToUse}
    matrix <Val00, Val01, Val02,
        Val10, Val11, Val12,
        Val20, Val21, Val22,
        Val30, Val31, Val32>
  }
}  
#end

 
 
#version testPart_Inc_Temp;
#end  