# Ray Tracing Algorithms

In this notebook, we will:

- Define a simple 3D scene (spheres, a plane, point lights, and a camera).
- Generate primary rays from a pinhole camera, one per pixel.
- Compute ray–object intersections analytically.
- Shade hit points with the Phong model, including shadows.
- Add recursive mirror reflections.
- Render the image and collect performance analytics (rays cast, hit rate, cost per recursion depth).
- Introduce a bounding-box test as the first step toward acceleration structures.

---

## 1. Setup

```python
import time
from dataclasses import dataclass, field

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams["figure.figsize"] = (8, 6)
plt.rcParams["figure.dpi"] = 120

RANDOM_SEED = 42
rng = np.random.default_rng(RANDOM_SEED)

EPS = 1e-4        # offset to avoid self-intersection ("shadow acne")
MAX_DEPTH = 3     # maximum number of recursive reflection bounces
```

---

## 2. Vector helpers

```python
def normalize(v):
    return v / np.linalg.norm(v)
```

---

## 3. Scene description

A ray tracer needs geometry (objects with an intersection test), materials (how a surface responds to light), and lights. We represent each as a small dataclass.

```python
@dataclass
class Material:
    color: np.ndarray
    ambient: float = 0.1
    diffuse: float = 0.7
    specular: float = 0.3
    shininess: float = 32.0
    reflectivity: float = 0.0   # 0 = matte, 1 = perfect mirror


@dataclass
class Sphere:
    center: np.ndarray
    radius: float
    material: Material

    def intersect(self, origin, direction):
        """Solve |origin + t*direction - center|^2 = radius^2 (direction is unit length)."""
        oc = origin - self.center
        b = np.dot(oc, direction)
        c = np.dot(oc, oc) - self.radius ** 2
        disc = b * b - c
        if disc < 0:
            return None
        sqrt_disc = np.sqrt(disc)
        for t in (-b - sqrt_disc, -b + sqrt_disc):
            if t > EPS:
                return t
        return None

    def normal(self, point):
        return normalize(point - self.center)


@dataclass
class Plane:
    point: np.ndarray
    normal_vec: np.ndarray
    material: Material

    def intersect(self, origin, direction):
        """Solve (origin + t*direction - point) . n = 0."""
        denom = np.dot(self.normal_vec, direction)
        if abs(denom) < 1e-8:          # ray parallel to plane
            return None
        t = np.dot(self.point - origin, self.normal_vec) / denom
        return t if t > EPS else None

    def normal(self, point):
        return self.normal_vec


@dataclass
class PointLight:
    position: np.ndarray
    intensity: float = 1.0


@dataclass
class Scene:
    objects: list
    lights: list
    background: np.ndarray = field(default_factory=lambda: np.array([0.05, 0.07, 0.12]))
```

Build a concrete scene:

```python
matte_red = Material(color=np.array([0.9, 0.2, 0.2]))
matte_blue = Material(color=np.array([0.2, 0.4, 0.9]))
chrome = Material(color=np.array([0.9, 0.9, 0.9]), ambient=0.05, diffuse=0.2,
                  specular=0.8, shininess=64.0, reflectivity=0.6)
floor = Material(color=np.array([0.8, 0.8, 0.6]), reflectivity=0.2)

scene = Scene(
    objects=[
        Sphere(np.array([-1.2, 1.0, 0.0]), 1.0, matte_red),
        Sphere(np.array([1.3, 0.6, 1.0]), 0.6, chrome),
        Sphere(np.array([0.2, 0.4, -1.5]), 0.4, matte_blue),
        Plane(np.array([0.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]), floor),
    ],
    lights=[
        PointLight(np.array([-3.0, 5.0, -4.0]), intensity=1.0),
        PointLight(np.array([4.0, 4.0, -2.0]), intensity=0.5),
    ],
)
```

---

## 4. Ray–scene intersection

The core query is: given a ray, which object does it hit first?

```python
def closest_hit(scene, origin, direction):
    best_t, best_obj = np.inf, None
    for obj in scene.objects:
        t = obj.intersect(origin, direction)
        if t is not None and t < best_t:
            best_t, best_obj = t, obj
    return best_t, best_obj
```

This brute-force loop is O(number of objects) per ray. Section 9 discusses how to reduce it.

---

## 5. Shadows

A point is in shadow if a ray from that point toward a light hits something before reaching the light.

```python
def in_shadow(scene, point, light_pos):
    to_light = light_pos - point
    dist = np.linalg.norm(to_light)
    direction = to_light / dist
    t, _ = closest_hit(scene, point + EPS * direction, direction)
    return t < dist
```

---

## 6. Shading (Phong model)

For each unshadowed light we add a diffuse term (Lambert's cosine law) and a specular term (highlight from the reflected light direction).

```python
def shade(scene, point, normal, view_dir, mat):
    color = mat.color * mat.ambient
    for light in scene.lights:
        if in_shadow(scene, point, light.position):
            continue
        L = normalize(light.position - point)          # direction to light
        diff = max(np.dot(normal, L), 0.0)
        R = 2.0 * np.dot(normal, L) * normal - L        # reflected light direction
        spec = max(np.dot(R, -view_dir), 0.0) ** mat.shininess
        color = color + light.intensity * (
            mat.diffuse * diff * mat.color + mat.specular * spec * np.ones(3)
        )
    return color
```

---

## 7. Recursive tracing with reflections

If a surface is reflective, we spawn a secondary ray in the mirror direction and blend its color in. Recursion stops at `MAX_DEPTH`.

```python
def new_stats():
    return {"rays": 0, "hits": 0, "rays_per_depth": [0] * (MAX_DEPTH + 1)}


def trace(scene, origin, direction, depth, stats):
    stats["rays"] += 1
    stats["rays_per_depth"][depth] += 1

    t, obj = closest_hit(scene, origin, direction)
    if obj is None:
        return scene.background

    stats["hits"] += 1
    point = origin + t * direction
    normal = obj.normal(point)
    color = shade(scene, point, normal, direction, obj.material)

    r = obj.material.reflectivity
    if r > 0.0 and depth < MAX_DEPTH:
        reflected = direction - 2.0 * np.dot(direction, normal) * normal
        bounce = trace(scene, point + EPS * normal, normalize(reflected),
                       depth + 1, stats)
        color = (1.0 - r) * color + r * bounce

    return color
```

---

## 8. Camera, rendering, and analytics

The camera generates one primary ray per pixel by mapping pixel coordinates onto an image plane located in front of a pinhole.

```python
def render(scene, width=320, height=240, fov_deg=45.0,
           cam_pos=np.array([0.0, 1.5, -6.0]),
           look_at=np.array([0.0, 0.5, 0.0])):
    forward = normalize(look_at - cam_pos)
    right = normalize(np.cross(forward, np.array([0.0, 1.0, 0.0])))
    up = np.cross(right, forward)

    aspect = width / height
    scale = np.tan(np.radians(fov_deg) / 2.0)

    img = np.zeros((height, width, 3))
    stats = new_stats()

    start = time.perf_counter()
    for j in range(height):
        for i in range(width):
            px = (2.0 * (i + 0.5) / width - 1.0) * aspect * scale
            py = (1.0 - 2.0 * (j + 0.5) / height) * scale
            d = normalize(forward + px * right + py * up)
            img[j, i] = trace(scene, cam_pos, d, 0, stats)
    elapsed = time.perf_counter() - start

    return np.clip(img, 0.0, 1.0), stats, elapsed
```

Render and display:

```python
W, H = 320, 240
image, stats, elapsed = render(scene, width=W, height=H)

plt.figure()
plt.imshow(image)
plt.axis("off")
plt.title("Ray-Traced Scene (Phong + Shadows + Reflections)")
plt.show()
```

Analytics:

```python
print(f"Render time:        {elapsed:.2f} s")
print(f"Pixels:             {W * H:,}")
print(f"Total rays cast:    {stats['rays']:,}")
print(f"Rays per pixel:     {stats['rays'] / (W * H):.2f}")
print(f"Hit rate:           {stats['hits'] / stats['rays']:.1%}")
print(f"Time per ray:       {elapsed / stats['rays'] * 1e6:.2f} µs")
print(f"Rays per depth:     {stats['rays_per_depth']}")

plt.figure()
plt.bar(range(MAX_DEPTH + 1), stats["rays_per_depth"])
plt.xlabel("Recursion depth")
plt.ylabel("Number of rays")
plt.title("Rays Cast per Recursion Depth")
plt.xticks(range(MAX_DEPTH + 1))
plt.show()
```

Expected pattern: depth 0 (primary rays) dominates, and each level of reflection casts fewer rays because only reflective surfaces spawn children.

---

## 9. Acceleration: bounding volumes (optional)

Testing every ray against every object is the main cost. A common first step is to wrap objects in axis-aligned bounding boxes (AABBs) and skip whole groups when a ray misses their box. Full acceleration structures (BVH, k-d trees, uniform grids) build on this idea.

The slab test below checks a ray against an AABB:

```python
def aabb_hit(bmin, bmax, origin, direction):
    inv = 1.0 / np.where(direction == 0.0, 1e-12, direction)
    t1 = (bmin - origin) * inv
    t2 = (bmax - origin) * inv
    t_near = np.max(np.minimum(t1, t2))
    t_far = np.min(np.maximum(t1, t2))
    return t_far >= max(t_near, 0.0)

# Quick check: a ray aimed at the red sphere should hit its bounding box
red_center = scene.objects[0].center
r = scene.objects[0].radius
ray_o = np.array([0.0, 1.0, -6.0])
ray_d = normalize(red_center - ray_o)
print("Hits bounding box:", aabb_hit(red_center - r, red_center + r, ray_o, ray_d))
```

To measure the benefit, count intersection tests with and without the box pre-check across a full render. This is a good extension exercise.

---

## 10. Summary and next steps

In this notebook we:

- Built a scene from spheres, a plane, and point lights.
- Generated primary rays from a pinhole camera.
- Computed analytic ray–sphere and ray–plane intersections.
- Added shadows and the Phong shading model.
- Implemented recursive reflections with a depth limit.
- Collected analytics on ray counts, hit rate, and cost per recursion depth.
- Introduced bounding-box tests as a first acceleration step.

Ideas to extend this analysis:

- Vectorize intersection tests with NumPy to replace the per-pixel Python loop.
- Add refraction using Snell's law and Fresnel equations for glass-like materials.
- Implement anti-aliasing by averaging multiple jittered samples per pixel.
- Build a BVH and compare intersection-test counts against brute force.
- Add area lights for soft shadows, or a Monte Carlo path tracer for global illumination.
- Compare render time and noise across resolutions and sample counts.
