# shadertoy球体渲染 + 简单phong光照模型

```javascript
const int maxMarchingSteps = 255;
const float minDistance = 0.0;
const float maxDistance = 100.0;
const float epsilon = 0.00001;

float sphereSDF(vec3 point){
    return length(point) - 1.0;
}

float sceneSDF(vec3 point){
    return sphereSDF(point);
}

float shortestDistanceToSurface(vec3 eye, vec3 direction, float start, float end){
    float depth = start;
    for (int i = 0; i < maxMarchingSteps; ++i){
        float distance = sceneSDF(eye + depth * direction);
        if (distance < epsilon){
            return depth;
        }
        depth += distance;
        if (depth > end){
            return end;
        }
    }
}

vec3 getDirection(float fieldOfView, vec2 size, vec2 fragCoord){
    vec2 xy = fragCoord - size / 2.0;
    float z = size.y / tan(radians(fieldOfView) / 2.0);
    return normalize(vec3(xy, -z));
}

vec3 estimateNormal(vec3 point){
    return normalize(vec3(
    sceneSDF(vec3(point.x + epsilon, point.y, point.z)) - sceneSDF(vec3(point.x - epsilon, point.y, point.z)),
    sceneSDF(vec3(point.x, point.y + epsilon, point.z)) - sceneSDF(vec3(point.x, point.y - epsilon, point.z)),
    sceneSDF(vec3(point.x, point.y, point.z + epsilon)) - sceneSDF(vec3(point.x, point.y, point.z - epsilon))
    ));
}

vec3 phongContributeForLight(vec3 diffuseColor, vec3 directColor, float alpha, vec3 pos, vec3 eye, vec3 lightPos, vec3 lightIntensity){
    vec3 normal = normalize(estimateNormal(pos));
    vec3 left = normalize(lightPos - pos);
    vec3 vertical = normalize(eye - pos);
    vec3 right = normalize(reflect(-left, normal));
    float dotLeft2Normal = dot(left, normal);
    float dotRight2Vertical = dot(right, vertical);
    if (dotLeft2Normal < 0.0){
        // 光源照射盲点
        return vec3(0.0, 0.0, 0.0);
    }
    if (dotRight2Vertical < 0.0){
        // 仅有漫反射
        return lightIntensity * (diffuseColor * dotLeft2Normal);
    }
    return lightIntensity * (diffuseColor * dotLeft2Normal + pow(dotRight2Vertical, alpha));
}

vec3 phongIllumination(vec3 environmentColor, vec3 diffuseColor, vec3 directColor, vec3 eye, float alpha, vec3 pos){
    const vec3 ambientLight = vec3(1.0, 1.0, 1.0) * 0.5;
    vec3 color = environmentColor * ambientLight;
    vec3 lightPos1= vec3(4.0 * sin(iTime), 2.0, 4.0 * cos(iTime));
    vec3 lightIntensity1 = vec3(0.4, 0.4, 0.4);
    color += phongContributeForLight(diffuseColor, directColor, alpha, pos, eye, lightPos1, lightIntensity1);
    vec3 lightPos2 = vec3(2.0 * sin(0.37 * iTime), 2.0 * cos(0.37 * iTime), 2.0);
    vec3 lightIntensity2 = vec3(0.4, 0.4, 0.4);
    color += phongContributeForLight(diffuseColor, directColor, alpha, pos, eye, lightPos2, lightIntensity2);
    return color;
}

void mainImage(out vec4 fragColor, in vec2 fragCoord){
    vec3 rayOrigin = vec3(0.0, 0.0, 5.0);
    vec3 rayDirection = getDirection(45.0, iResolution.xy, fragCoord);
    float distance = shortestDistanceToSurface(rayOrigin, rayDirection, minDistance, maxDistance);
    if (distance > maxDistance - epsilon){
        fragColor = vec4(0.0, 0.0, 0.0, 0.0);
        return;
    }
    vec3 pos = rayOrigin + distance * rayDirection;
    vec3 environmentColor = vec3(0.2, 0.2, 0.2);
    vec3 diffuseColor = vec3(0.7, 0.2, 0.2);
    vec3 directColor = vec3(1.0, 1.0, 1.0);
    float shininess = 10.0;
    vec3 color = phongIllumination(environmentColor, diffuseColor, directColor, rayOrigin, shininess, pos);
    fragColor = vec4(color, 1.0);
    return;
}
```

# shadertoy高斯模糊

```javascript
float normpdf(in float x, in float sigma)
{
	return 0.39894*exp(-0.5*x*x/(sigma*sigma))/sigma;
}

void mainImage(out vec4 fragColor, in vec2 fragCoord){
    vec3 c = texture(iChannel0, fragCoord.xy / iResolution.xy).rgb;
    const int count = 11;
    const int middle = (count - 1) / 2;
    float kernel[count];
    vec3 color = vec3(0.0, 0.0, 0.0);
    float sigma = 5.0;
    float z = 0.0;
    for (int i = 0; i <= middle; i++){
        kernel[middle - i] = kernel[middle + i] = normpdf(float(i), sigma);
    }
    for (int i = 0; i < count; i++){
        z += kernel[i];
    }
    for (int i = -middle; i <= middle; i++){
        for (int j = -middle; j <= middle; j++){
            color += kernel[middle + j] * kernel[middle + i] * texture(iChannel0, (fragCoord.xy + vec2(float(i), float(j))) / iResolution.xy).rgb;
        }
    }
    fragColor = vec4(color/(z * z), 1.0);
}
```

