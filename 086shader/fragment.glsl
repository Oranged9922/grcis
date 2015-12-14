//---------------------------------------------------------------
// Simple fragment shader utilizing single texture
//   normal vector, single point light-source
//   @version $Rev$
//---------------------------------------------------------------
#version 330

in vec2 varTexCoords;
in vec3 varNormal;
in vec3 varWorld;
in vec3 varColor;

uniform vec3 globalAmbient;
uniform vec3 lightColor;
uniform vec3 lightPosition;
uniform vec3 eyePosition;
uniform vec3 Ka;
uniform vec3 Kd;
uniform vec3 Ks;
uniform float shininess;
uniform bool useTexture;
uniform bool globalColor;
uniform sampler2D texWater;

out vec4 fragColor;

void main ()
{
  vec3 P = varWorld;
  vec3 N = normalize( varNormal );
  vec3 L = normalize( lightPosition - P );
  vec3 V = normalize( eyePosition - P );
  vec3 H = normalize( L + V );

  float cosb = pow( max( dot( N, H ), 0.0 ), shininess );
  float cosa = dot( N, L );
  if ( cosa < 0.0 ) cosb = 0.0;
  cosa = max( cosa, 0.0 );

  vec3 ka, kd;
  if ( useTexture )
  {
    ka = kd = vec3( texture2D( texWater, varTexCoords ) );
  }
  else
    if ( globalColor )
    {
      ka = Ka;
      kd = Kd;
    }
    else
    {
      ka = kd = varColor;
    }
  vec3 ambient  = ka * globalAmbient;
  vec3 diffuse  = kd * lightColor * cosa;
  vec3 specular = Ks * lightColor * cosb;

  fragColor.rgb = ambient + diffuse + specular;
  fragColor.a   = 1.0;
}
