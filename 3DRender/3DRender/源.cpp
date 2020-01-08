#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<assert.h>
#include<windows.h>
#include<tchar.h>
#include<iostream>
#include <vector>
using namespace std;
typedef unsigned int UINT;

// 数学函数

typedef struct {
	float matrix[4][4];
} matrix_t;

typedef struct {
	float x, y, z, w;
} vector_t, point_t;

// 插值 t[0, 1]
float interp(float x1, float x2, float t) {
	return x1 + (x2 - x1) * t;
}

int CMID(int x, int min, int max) { return (x < min) ? min : ((x > max) ? max : x); }


// 向量的模
float vector_length(const vector_t* v) {
	float sq = v->x * v->x + v->y * v->y + v->z * v->z;
	return (float)sqrt(sq);
}

// 向量点乘
float vector_dotproduct(const vector_t* x, const vector_t* y) {
	return x->x * y->x + x->y * y->y + x->z * y->z;
}

// 向量叉乘
void vector_crossproduct(vector_t* z, const vector_t* x, const vector_t* y) {
	float x_result, y_result, z_result;
	x_result = x->y * y->z - x->z * y->y;
	y_result = x->z * y->x - x->x * y->z;
	z_result = x->x * y->y - x->y * y->x;
	z->x = x_result;
	z->y = y_result;
	z->z = z_result;
	z->w = 1.0f;
}

// 向量加法 vector = vector1 + vector2
void vector_add(vector_t* vector, const vector_t* vector1, const vector_t* vector2) {
	vector->x = vector1->x + vector2->x;
	vector->y = vector1->y + vector2->y;
	vector->z = vector1->z + vector2->y;
}

// 向量减法 vector = vector1 - vector2
void vector_sub(vector_t* vector, const vector_t* vector1, const vector_t* vector2) {
	vector->x = vector1->x - vector2->x;
	vector->y = vector1->y - vector2->y;
	vector->z = vector1->z - vector2->y;
}

// 向量的插值
void vector_interp(vector_t* z, const vector_t* x1, const vector_t* x2, float t) {
	z->x = interp(x1->x, x2->x, t);
	z->y = interp(x1->y, x2->y, t);
	z->z = interp(x1->z, x2->z, t);
	z->w = 1.0f;
}

// 向量归一化
void vector_normalize(vector_t* v) {
	float length = vector_length(v);
	if (length != 0.0f) {
		float inv = 1.0f / length;
		v->x *= inv;
		v->y *= inv;
		v->z *= inv;
	}
}



// 矩阵相加 c = a + b
void matrix_add(matrix_t* c, const matrix_t* a, const matrix_t* b) {
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			c->matrix[i][j] = a->matrix[i][j] + b->matrix[i][j];
		}
	}
}

// 矩阵相减 c = a - b
void matrix_sub(matrix_t* c, const matrix_t* a, const matrix_t* b) {
	int i, j;
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++)
			c->matrix[i][j] = a->matrix[i][j] - b->matrix[i][j];
	}
}

// 矩阵相乘 c = a * b
void matrix_mul(matrix_t* c, const matrix_t* a, const matrix_t* b) {
	matrix_t z;
	int i, j;
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			z.matrix[j][i] = (a->matrix[j][0] * b->matrix[0][i]) +
				(a->matrix[j][1] * b->matrix[1][i]) +
				(a->matrix[j][2] * b->matrix[2][i]) +
				(a->matrix[j][3] * b->matrix[3][i]);
		}
	}
	c[0] = z;
}

// 矩阵缩放 c = a * f
void matrix_scale(matrix_t* c, const matrix_t* a, float f) {
	int i, j;
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++)
			c->matrix[i][j] = a->matrix[i][j] * f;
	}
}

// 向量右乘矩阵，得到一个新的向量 y = x * m
void matrix_apply(vector_t* y, const vector_t* x, const matrix_t* m) {
	float X = x->x, Y = x->y, Z = x->z, W = x->w;
	y->x = X * m->matrix[0][0] + Y * m->matrix[1][0] + Z * m->matrix[2][0] + W * m->matrix[3][0];
	y->y = X * m->matrix[0][1] + Y * m->matrix[1][1] + Z * m->matrix[2][1] + W * m->matrix[3][1];
	y->z = X * m->matrix[0][2] + Y * m->matrix[1][2] + Z * m->matrix[2][2] + W * m->matrix[3][2];
	y->w = X * m->matrix[0][3] + Y * m->matrix[1][3] + Z * m->matrix[2][3] + W * m->matrix[3][3];
}

// 将矩阵设置为零矩阵
void matrix_set_zero(matrix_t* m) {
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			m->matrix[i][j] = 0.0f;
		}
	}
}

// 将矩阵设置为单位矩阵
void matrix_set_identity(matrix_t* m) {
	matrix_set_zero(m);
	m->matrix[0][0] = m->matrix[1][1] = m->matrix[2][2] = m->matrix[3][3] = 1.0f;
}

// 变换：平移矩阵
void matrix_set_translate(matrix_t* m, float x, float y, float z) {
	matrix_set_identity(m);
	m->matrix[3][0] = x;
	m->matrix[3][1] = y;
	m->matrix[3][2] = z;
}

// 变换：缩放矩阵
void matrix_set_scale(matrix_t* m, float x, float y, float z) {
	matrix_set_identity(m);
	m->matrix[0][0] = x;
	m->matrix[1][1] = y;
	m->matrix[2][2] = z;
}

// 变换：旋转矩阵
void matrix_set_rotate(matrix_t* m, float x, float y, float z, float theta) {
	float qsin = (float)sin(theta * 0.5);
	float qcos = (float)cos(theta * 0.5);
	vector_t vec = { x, y, z, 1.0f };
	float w = qcos;
	vector_normalize(&vec);
	x = vec.x * qsin;
	y = vec.y * qsin;
	z = vec.z * qsin;
	m->matrix[0][0] = 1 - 2 * y * y - 2 * z * z;
	m->matrix[1][0] = 2 * x * y - 2 * w * z;
	m->matrix[2][0] = 2 * x * z + 2 * w * y;
	m->matrix[0][1] = 2 * x * y + 2 * w * z;
	m->matrix[1][1] = 1 - 2 * x * x - 2 * z * z;
	m->matrix[2][1] = 2 * y * z - 2 * w * x;
	m->matrix[0][2] = 2 * x * z - 2 * w * y;
	m->matrix[1][2] = 2 * y * z + 2 * w * x;
	m->matrix[2][2] = 1 - 2 * x * x - 2 * y * y;
	m->matrix[0][3] = m->matrix[1][3] = m->matrix[2][3] = 0.0f;
	m->matrix[3][0] = m->matrix[3][1] = m->matrix[3][2] = 0.0f;
	m->matrix[3][3] = 1.0f;
}

// 设置摄像机（得到视图变换矩阵，即仿射变换）
void matrix_set_lookat(matrix_t* m, const vector_t* eye, const vector_t* at, const vector_t* up) {
	vector_t xaxis, yaxis, zaxis;

	// 计算摄像机本地坐标系
	vector_sub(&zaxis, at, eye);
	vector_normalize(&zaxis);
	vector_crossproduct(&xaxis, up, &zaxis);
	vector_normalize(&xaxis);
	vector_crossproduct(&yaxis, &zaxis, &xaxis);
	vector_normalize(&yaxis);

	// 计算视图变换矩阵
	m->matrix[0][0] = xaxis.x;
	m->matrix[1][0] = xaxis.y;
	m->matrix[2][0] = xaxis.z;
	m->matrix[3][0] = -vector_dotproduct(&xaxis, eye);

	m->matrix[0][1] = yaxis.x;
	m->matrix[1][1] = yaxis.y;
	m->matrix[2][1] = yaxis.z;
	m->matrix[3][1] = -vector_dotproduct(&yaxis, eye);

	m->matrix[0][2] = zaxis.x;
	m->matrix[1][2] = zaxis.y;
	m->matrix[2][2] = zaxis.z;
	m->matrix[3][2] = -vector_dotproduct(&zaxis, eye);

	m->matrix[0][3] = m->matrix[1][3] = m->matrix[2][3] = 0.0f;
	m->matrix[3][3] = 1.0f;
}

// 根据屏幕宽高比和镜头视角大小设置透视投影矩阵，还有一种方案是根据屏幕边界计算
void matrix_set_perspective(matrix_t* m, float fovy, float aspect, float zn, float zf) {
	float fax = 1.0f / (float)tan(fovy * 0.5f);
	matrix_set_zero(m);
	m->matrix[0][0] = (float)(fax / aspect);
	m->matrix[1][1] = (float)(fax);
	m->matrix[2][2] = zf / (zf - zn);
	m->matrix[3][2] = -zn * zf / (zf - zn);
	m->matrix[2][3] = 1;
}

// 结构体内的矩阵都是变换矩阵，transform 是最终的变换矩阵
typedef struct {
	matrix_t world;
	matrix_t view;
	matrix_t projection;
	matrix_t transform; // transform = world * view * projection
	float screen_width, screen_height;
} transform_t;

// 变换矩阵刷新
void transform_update(transform_t* transform) {
	matrix_t m;
	matrix_mul(&m, &(transform->world), &(transform->view));
	matrix_mul(&(transform->transform), &m, &(transform->projection));
}

// 初始化，设置屏幕长宽
void transform_init(transform_t* transform, int width, int height) {
	float aspect = (float)width / ((float)height);
	matrix_set_identity(&transform->world);
	matrix_set_identity(&transform->view);
	matrix_set_perspective(&transform->projection, 3.1415926f * 0.5f, aspect, 1.0f, 500.0f);
	transform->screen_width = (float)width;
	transform->screen_height = (float)height;
	transform_update(transform);
}

// 将世界坐标系的向量变换到以摄像机为中心的坐标系中
void transform_apply(vector_t* vector_in_camera, const vector_t* vector_in_world, const transform_t* ts) {
	matrix_apply(vector_in_camera, vector_in_world, &ts->transform);
}

// 归一化，得到屏幕坐标（屏幕上的像素坐标）
void transform_homogenize(const transform_t* ts, vector_t* y, const vector_t* x) {
	float rhw = 1.0f / x->w;
	y->x = (x->x * rhw + 1.0f) * ts->screen_width * 0.5f;
	y->y = (1.0f - x->y * rhw) * ts->screen_height * 0.5f;
	y->z = x->z * rhw;
	y->w = 1.0f;
}

// 绘制图元所需的几种数据
typedef struct {
	float r, g, b;
} color_t;

typedef struct {
	point_t pos;
	color_t color;
	float rhw;
} vertex_t;

typedef struct {
	vertex_t v, v1, v2;
} edge_t;

typedef struct {
	float top, bottom;
	edge_t left, right;
} trapezoid_t;

typedef struct {
	vertex_t v, step;
	int x, y, scan_width;
} scanline_t;

void vertex_rhw_init(vertex_t* v) {
	float rhw = 1.0f / v->pos.w;
	v->rhw = rhw;
	v->color.r *= rhw;
	v->color.g *= rhw;
	v->color.b *= rhw;
}

// 顶点的插值计算
void vertex_interp(vertex_t* y, const vertex_t* x1, const vertex_t* x2, float t) {
	vector_interp(&y->pos, &x1->pos, &x2->pos, t);
	y->color.r = interp(x1->color.r, x2->color.r, t);
	y->color.g = interp(x1->color.g, x2->color.g, t);
	y->color.b = interp(x1->color.b, x2->color.b, t);
	y->rhw = interp(x1->rhw, x2->rhw, t);
}

void vertex_division(vertex_t* y, const vertex_t* x1, const vertex_t* x2, float width) {
	float inv = 1.0f / width;
	y->pos.x = (x2->pos.x - x1->pos.x) * inv;
	y->pos.y = (x2->pos.y - x1->pos.y) * inv;
	y->pos.z = (x2->pos.z - x1->pos.z) * inv;
	y->pos.w = (x2->pos.w - x1->pos.w) * inv;
	y->color.r = (x2->color.r - x1->color.r) * inv;
	y->color.g = (x2->color.g - x1->color.g) * inv;
	y->color.b = (x2->color.b - x1->color.b) * inv;
	y->rhw = (x2->rhw - x1->rhw) * inv;
}


// 顶点的相加
void vertex_add(vertex_t* y, const vertex_t* x) {
	y->pos.x += x->pos.x;
	y->pos.y += x->pos.y;
	y->pos.z += x->pos.z;
	y->pos.w += x->pos.w;
	y->rhw += x->rhw;
	y->color.r += x->color.r;
	y->color.g += x->color.g;
	y->color.b += x->color.b;
}

int trapezoid_init_triangle(trapezoid_t* trap, const vertex_t* v1, const vertex_t* v2, const vertex_t* v3) {
	const vertex_t* v;
	if (v1->pos.y > v2->pos.y) {
		v = v1;
		v1 = v2;
		v2 = v;
	}
	if (v1->pos.y > v3->pos.y) {
		v = v1;
		v1 = v3;
		v3 = v;
	}
	if (v2->pos.y > v3->pos.y) {
		v = v2;
		v2 = v3;
		v3 = v;
	}
	if (v1->pos.y == v2->pos.y && v1->pos.y == v3->pos.y) {
		return 0;
	}
	if (v1->pos.x == v2->pos.x && v1->pos.x == v3->pos.x) {
		return 0;
	}
	if (v1->pos.y == v2->pos.y) {
		if (v1->pos.x > v2->pos.x) {
			v = v1;
			v1 = v2;
			v2 = v;
		}
		trap[0].top = v1->pos.y;
		trap[0].bottom = v3->pos.y;
		trap[0].left.v1 = *v1;
		trap[0].left.v2 = *v3;
		trap[0].right.v1 = *v2;
		trap[0].right.v2 = *v3;
		return trap[0].top < trap[0].bottom ? 1 : 0;
	}
	if (v2->pos.y == v3->pos.y) {
		if (v2->pos.x > v3->pos.x) {
			v = v2;
			v2 = v3;
			v3 = v;
		}
		trap[0].top = v1->pos.y;
		trap[0].bottom = v2->pos.y;
		trap[0].left.v1 = *v1;
		trap[0].left.v2 = *v2;
		trap[0].right.v1 = *v1;
		trap[0].right.v2 = *v3;
		return trap[0].top < trap[0].bottom ? 1 : 0;
	}
	trap[0].top = v1->pos.y;
	trap[0].bottom = v2->pos.y;
	trap[1].top = v2->pos.y;
	trap[1].bottom = v3->pos.y;

	float k = (v3->pos.y - v1->pos.y) / (v2->pos.y - v1->pos.y);
	float x = v1->pos.x + (v2->pos.x - v1->pos.x) * k;

	if (x <= v3->pos.x) {
		trap[0].left.v1 = *v1;
		trap[0].left.v2 = *v2;
		trap[0].right.v1 = *v1;
		trap[0].right.v2 = *v3;
		trap[1].left.v1 = *v2;
		trap[1].left.v2 = *v3;
		trap[1].right.v1 = *v1;
		trap[1].right.v2 = *v3;

	}
	else {
		trap[0].left.v1 = *v1;
		trap[0].left.v2 = *v3;
		trap[0].right.v1 = *v1;
		trap[0].right.v2 = *v2;
		trap[1].left.v1 = *v1;
		trap[1].left.v2 = *v3;
		trap[1].right.v1 = *v2;
		trap[1].right.v2 = *v3;
	}
	return 2;
}

// 按参数b值计算出直线 y=m 与四边形的交点，存储在四边形结构体中
void trapezoid_edge_interp(trapezoid_t* trap, float m) {
	// 处理左边的线段
	float s1 = trap->left.v2.pos.y - trap->left.v1.pos.y;
	float t1 = (m - trap->left.v1.pos.y) / s1;
	vertex_interp(&(trap->left.v), &(trap->left.v1), &(trap->left.v2), t1);

	// 处理右边的线段
	float s2 = trap->right.v2.pos.y - trap->right.v1.pos.y;
	float t2 = (m - trap->right.v1.pos.y) / s2;
	vertex_interp(&(trap->right.v), &(trap->right.v1), &(trap->right.v2), t2);
}

// 根据上一步计算出的交点初始化扫描线和步长
void trapezoid_init_scanline(trapezoid_t* trap, scanline_t* scanline, int y) {
	// 初始化扫描线
	float width = trap->right.v.pos.x - trap->left.v.pos.x;
	scanline->x = (int)(trap->left.v.pos.x + 0.5f);
	scanline->scan_width = (int)(trap->right.v.pos.x + 0.5f) - scanline->x;
	scanline->y = y;
	scanline->v = trap->left.v;
	if (trap->left.v.pos.x >= trap->right.v.pos.x) {
		scanline->scan_width = 0;
	}
	// 计算步长
	vertex_division(&scanline->step, &trap->left.v, &trap->right.v, width);
}

// 渲染设备
typedef struct {
	transform_t transform;      // 坐标变换器
	int width;                  // 窗口宽度
	int height;                 // 窗口高度
	UINT** framebuffer;      // 像素缓存：framebuffer[y] 代表第 y行
	float** zbuffer;            // 深度缓存：zbuffer[y] 为第 y行指针
	float max_u;                // 纹理最大宽度：tex_width - 1
	float max_v;                // 纹理最大高度：tex_height - 1
	int render_state;           // 渲染状态
	UINT background;         // 背景颜色
	UINT foreground;         // 线框颜色
}	device_t;


#define RENDER_STATE_WIREFRAME      1		// 渲染线框
#define RENDER_STATE_COLOR          2		// 渲染颜色

// 设备初始化，fb为外部帧缓存，非 NULL 将引用外部帧缓存（每行 4字节对齐）
void device_init(device_t * device, int width, int height, void* fb) {
	int need = sizeof(void*) * (height * 2 + 1024) + width * height * 8;
	char* ptr = (char*)malloc(need + 64);
	char* framebuf, * zbuf;
	int j;
	assert(ptr);
	device->framebuffer = (UINT * *)ptr;
	device->zbuffer = (float**)(ptr + sizeof(void*) * height);
	ptr += sizeof(void*) * height * 2;
	ptr += sizeof(void*) * 1024;
	framebuf = (char*)ptr;
	zbuf = (char*)ptr + width * height * 4;
	ptr += width * height * 8;
	if (fb != NULL) framebuf = (char*)fb;
	for (j = 0; j < height; j++) {
		device->framebuffer[j] = (UINT*)(framebuf + width * 4 * j);
		device->zbuffer[j] = (float*)(zbuf + width * 4 * j);
	}

	device->max_u = 1.0f;
	device->max_v = 1.0f;
	device->width = width;
	device->height = height;
	device->background = 0xc0c0c0;
	device->foreground = 0;
	transform_init(&device->transform, width, height);
	device->render_state = RENDER_STATE_WIREFRAME;
}

// 删除设备
void device_destroy(device_t* device) {
	if (device->framebuffer)
		free(device->framebuffer);
	device->framebuffer = NULL;
	device->zbuffer = NULL;
}

// 设置当前纹理
void device_set_texture(device_t* device, void* bits, long pitch, int w, int h) {
	char* ptr = (char*)bits;
	int j;
	assert(w <= 1024 && h <= 1024);
	device->max_u = (float)(w - 1);
	device->max_v = (float)(h - 1);
}

// 清空 framebuffer 和 zbuffer
void device_clear(device_t* device, int mode) {
	int y, x, height = device->height;
	for (y = 0; y < device->height; y++) {
		UINT* dst = device->framebuffer[y];
		UINT cc = (height - 1 - y) * 230 / (height - 1);
		cc = (cc << 16) | (cc << 8) | cc;
		if (mode == 0) cc = device->background;
		for (x = device->width; x > 0; dst++, x--) dst[0] = cc;
	}
	for (y = 0; y < device->height; y++) {
		float* dst = device->zbuffer[y];
		for (x = device->width; x > 0; dst++, x--) dst[0] = 0.0f;
	}
}

// 画点
void device_pixel(device_t* device, int x, int y, UINT32 color) {
	if (((UINT)x) < (UINT)device->width && ((UINT)y) < (UINT)device->height) {
		device->framebuffer[y][x] = color;
	}
}

// 画直线
void device_draw_line(device_t* device, int x1, int y1, int x2, int y2, UINT color) {

	int x, y;

	if (x1 == x2 && y1 == y2) {
		device_pixel(device, x1, y1, color);
	}
	else if (x1 == x2) {
		int sign = y1 <= y2 ? 1 : -1;
		for (y = y1; y != y2; y += sign) {
			device_pixel(device, x1, y, color);
		}
		device_pixel(device, x2, y2, color);
	}
	else if (y1 == y2) {
		int sign = x1 <= x2 ? 1 : -1;
		for (x = x1; x != x2; x += sign) {
			device_pixel(device, x, y1, color);
		}
		device_pixel(device, x2, y2, color);
	}
	else {
		// 以上是三种斜率不存在的情况，以下是直线斜率存在的情况
		int dx = x1 < x2 ? x2 - x1 : x1 - x2;
		int dy = y1 < y2 ? y2 - y1 : y1 - y2;
		if (dy <= dx) {
			int rem = 0;
			// 直线斜率的绝对值小于等于1
			if (x1 > x2) {
				// 将 x 较小者设置为起始点
				x = x1;
				x1 = x2;
				x2 = x;
				y = y1;
				y1 = y2;
				y2 = y;
			}
			for (x = x1, y = y1; x <= x2; ++x) {
				device_pixel(device, x, y, color);
				rem += dy;
				if (rem >= dx) {
					rem -= dx;
					y += (y1 <= y2) ? 1 : -1;
					device_pixel(device, x, y, color);
				}
			}
			device_pixel(device, x2, y2, color);
		} else {
			int rem = 0;
			// 直线斜率的绝对值大于1，交换 x 和 y
			if (y1 > y2) {
				x = x1;
				x1 = x2;
				x2 = x;
				y = y1;
				y1 = y2;
				y2 = y;
			}
			for (x = x1, y = y1; y <= y2; ++y) {
				device_pixel(device, x, y, color);
				rem += dx;
				if (rem >= dy) {
					rem -= dy;
					x += (x1 <= x2 )? 1 : -1;
					device_pixel(device, x, y, color);
				}
			}
			device_pixel(device, x2, y2, color);
		}
	}
}

// 绘制扫描线
void device_draw_scanline(device_t* device, scanline_t* scanline) {
	UINT* framebuffer = device->framebuffer[scanline->y];
	float* zbuffer = device->zbuffer[scanline->y];
	int scan_x = scanline->x;
	int screen_width = device->width;
	int scan_width = scanline->scan_width;
	int render_state = device->render_state;
	for (; scan_width > 0; --scan_width, ++scan_x) {
		if (scan_x > 0 && scan_x < screen_width) {
			float rhw = scanline->v.rhw;
			// rhw 越大则 w 值越小，意味着越靠近摄像机，越优先被渲染
			if (rhw >= zbuffer[scan_x]) {
				float w = 1.0f / rhw;
				zbuffer[scan_x] = rhw;
				if (render_state && RENDER_STATE_COLOR) {
					float r = scanline->v.color.r * w;
					float g = scanline->v.color.g * w;
					float b = scanline->v.color.b * w;
					int R = (int)(r * 255.0f);
					int G = (int)(g * 255.0f);
					int B = (int)(b * 255.0f);
					R = CMID(R, 0, 255);
					G = CMID(G, 0, 255);
					B = CMID(B, 0, 255);
					framebuffer[scan_x] = (R << 16) | (G << 8) | (B);
				}
			}
		}
		// 迭代扫描的顶点
		vertex_add(&scanline->v, &scanline->step);
		if (scan_x >= screen_width) {
			break;
		}
	}
}

// 渲染四边形（由三角形图元分割而来）
void device_render_trap(device_t* device, trapezoid_t* trap) {
	scanline_t scanline;
	int top = (int)(trap->top + 0.5f);
	int bottom = (int)(trap->bottom + 0.5f);
	for (int i = top; i < bottom; ++i) {
		if (i >= 0 && i < device->height) {
			// 在屏幕范围内
			trapezoid_edge_interp(trap, (float)i + 0.5f); // 计算出四边形与水平渲染线的交点

			// 计算出交点则已知哪些顶点是在图元内的，防止污染图元外部的顶点
			trapezoid_init_scanline(trap, &scanline, i);
			device_draw_scanline(device, &scanline);
		}
		if (i >= device->height) {
			break;
		}
	}
}

// 检查齐次坐标同 cvv 的边界用于视锥裁剪
int transform_check_cvv(const vector_t* v) {
	float w = v->w;
	int check = 0;
	if (v->z < 0.0f) check |= 1;
	if (v->z > w) check |= 2;
	if (v->x < -w) check |= 4;
	if (v->x > w) check |= 8;
	if (v->y < -w) check |= 16;
	if (v->y > w) check |= 32;
	return check;
}

// 绘制图元（三角形）
void device_draw_primitive(device_t* device, const vertex_t* v1, const vertex_t* v2, const vertex_t* v3) {
	int render_state = device->render_state;
	point_t screen_pos1, screen_pos2, screen_pos3, view_pos1, view_pos2, view_pos3;

	// 先将给的三个顶点变换成视图坐标，存储在view_pos*中
	transform_apply(&view_pos1, &v1->pos, &device->transform);
	transform_apply(&view_pos2, &v2->pos, &device->transform);
	transform_apply(&view_pos3, &v3->pos, &device->transform);

	// 检查变换后的顶点坐标是否在相机范围内
	if (transform_check_cvv(&view_pos1) != 0) return;
	if (transform_check_cvv(&view_pos2) != 0) return;
	if (transform_check_cvv(&view_pos3) != 0) return;

	// 得到屏幕坐标（以像素计），存储在screen_pos*中
	transform_homogenize(&device->transform, &screen_pos1, &view_pos1);
	transform_homogenize(&device->transform, &screen_pos2, &view_pos2);
	transform_homogenize(&device->transform, &screen_pos3, &view_pos3);

	// 彩色绘制
	if (render_state & RENDER_STATE_COLOR) {
		vertex_t t1 = *v1, t2 = *v2, t3 = *v3;
		trapezoid_t traps[2];
		int n;

		t1.pos = screen_pos1;
		t2.pos = screen_pos2;
		t3.pos = screen_pos3;
		t1.pos.w = view_pos1.w;
		t2.pos.w = view_pos2.w;
		t3.pos.w = view_pos3.w;

		vertex_rhw_init(&t1);	// 初始化 w
		vertex_rhw_init(&t2);	// 初始化 w
		vertex_rhw_init(&t3);	// 初始化 w

		// 拆分三角形为0-2个梯形，并且返回可用梯形数量
		n = trapezoid_init_triangle(traps, &t1, &t2, &t3);
		// 渲染梯形
		if (n >= 1) device_render_trap(device, &traps[0]);
		if (n >= 2) device_render_trap(device, &traps[1]);
	}
	if (render_state & RENDER_STATE_WIREFRAME) {
		// 线框绘制
		device_draw_line(device, (int)screen_pos1.x, (int)screen_pos1.y, (int)screen_pos2.x, (int)screen_pos2.y, device->foreground);
		device_draw_line(device, (int)screen_pos1.x, (int)screen_pos1.y, (int)screen_pos3.x, (int)screen_pos3.y, device->foreground);
		device_draw_line(device, (int)screen_pos3.x, (int)screen_pos3.y, (int)screen_pos2.x, (int)screen_pos2.y, device->foreground);
	}
}

int screen_w, screen_h, screen_exit = 0;
int screen_mx = 0, screen_my = 0, screen_mb = 0;
int screen_keys[512];	// 当前键盘按下状态
static HWND screen_handle = NULL;		// 主窗口 HWND
static HDC screen_dc = NULL;			// 配套的 HDC
static HBITMAP screen_hb = NULL;		// DIB
static HBITMAP screen_ob = NULL;		// 老的 BITMAP
unsigned char* screen_fb = NULL;		// frame buffer
long screen_pitch = 0;

int screen_init(int w, int h, const TCHAR* title);	// 屏幕初始化
int screen_close(void);								// 关闭屏幕
void screen_dispatch(void);							// 处理消息
void screen_update(void);							// 显示 FrameBuffer

// win32 event handler
static LRESULT screen_events(HWND, UINT, WPARAM, LPARAM);

#ifdef _MSC_VER
#pragma comment(lib, "gdi32.lib")
#pragma comment(lib, "user32.lib")
#endif

// 初始化窗口并设置标题
int screen_init(int w, int h, const TCHAR* title) {
	WNDCLASS wc = { CS_BYTEALIGNCLIENT, (WNDPROC)screen_events, 0, 0, 0,
		NULL, NULL, NULL, NULL, _T("SCREEN3.1415926") };
	BITMAPINFO bi = { { sizeof(BITMAPINFOHEADER), w, -h, 1, 32, BI_RGB,
		w * h * 4, 0, 0, 0, 0 } };
	RECT rect = { 0, 0, w, h };
	int wx, wy, sx, sy;
	LPVOID ptr;
	HDC hDC;

	screen_close();

	wc.hbrBackground = (HBRUSH)GetStockObject(BLACK_BRUSH);
	wc.hInstance = GetModuleHandle(NULL);
	wc.hCursor = LoadCursor(NULL, IDC_ARROW);
	if (!RegisterClass(&wc)) return -1;

	screen_handle = CreateWindow(_T("SCREEN3.1415926"), title,
		WS_OVERLAPPED | WS_CAPTION | WS_SYSMENU | WS_MINIMIZEBOX,
		0, 0, 0, 0, NULL, NULL, wc.hInstance, NULL);
	if (screen_handle == NULL) return -2;

	screen_exit = 0;
	hDC = GetDC(screen_handle);
	screen_dc = CreateCompatibleDC(hDC);
	ReleaseDC(screen_handle, hDC);

	screen_hb = CreateDIBSection(screen_dc, &bi, DIB_RGB_COLORS, &ptr, 0, 0);
	if (screen_hb == NULL) return -3;

	screen_ob = (HBITMAP)SelectObject(screen_dc, screen_hb);
	screen_fb = (unsigned char*)ptr;
	screen_w = w;
	screen_h = h;
	screen_pitch = w * 4;

	AdjustWindowRect(&rect, GetWindowLong(screen_handle, GWL_STYLE), 0);
	wx = rect.right - rect.left;
	wy = rect.bottom - rect.top;
	sx = (GetSystemMetrics(SM_CXSCREEN) - wx) / 2;
	sy = (GetSystemMetrics(SM_CYSCREEN) - wy) / 2;
	if (sy < 0) sy = 0;
	SetWindowPos(screen_handle, NULL, sx, sy, wx, wy, (SWP_NOCOPYBITS | SWP_NOZORDER | SWP_SHOWWINDOW));
	SetForegroundWindow(screen_handle);

	ShowWindow(screen_handle, SW_NORMAL);
	screen_dispatch();

	memset(screen_keys, 0, sizeof(int) * 512);
	memset(screen_fb, 0, w * h * 4);

	return 0;
}

int screen_close(void) {
	if (screen_dc) {
		if (screen_ob) {
			SelectObject(screen_dc, screen_ob);
			screen_ob = NULL;
		}
		DeleteDC(screen_dc);
		screen_dc = NULL;
	}
	if (screen_hb) {
		DeleteObject(screen_hb);
		screen_hb = NULL;
	}
	if (screen_handle) {
		CloseWindow(screen_handle);
		screen_handle = NULL;
	}
	return 0;
}

static LRESULT screen_events(HWND hWnd, UINT msg,
	WPARAM wParam, LPARAM lParam) {
	switch (msg) {
	case WM_CLOSE: screen_exit = 1; break;
	case WM_KEYDOWN: screen_keys[wParam & 511] = 1; break;
	case WM_KEYUP: screen_keys[wParam & 511] = 0; break;
	default: return DefWindowProc(hWnd, msg, wParam, lParam);
	}
	return 0;
}

void screen_dispatch(void) {
	MSG msg;
	while (1) {
		if (!PeekMessage(&msg, NULL, 0, 0, PM_NOREMOVE)) break;
		if (!GetMessage(&msg, NULL, 0, 0)) break;
		DispatchMessage(&msg);
	}
}

void screen_update(void) {
	HDC hDC = GetDC(screen_handle);
	BitBlt(hDC, 0, 0, screen_w, screen_h, screen_dc, 0, 0, SRCCOPY);
	ReleaseDC(screen_handle, hDC);
	screen_dispatch();
}

// cvv 的 8 个顶点
vertex_t mesh[8] = {
	{ { -1, -1,  1, 1 }, { 1.0f, 0.2f, 0.2f }, 1 },
	{ {  1, -1,  1, 1 }, { 0.2f, 1.0f, 0.2f }, 1 },
	{ {  1,  1,  1, 1 }, { 0.2f, 0.2f, 1.0f }, 1 },
	{ { -1,  1,  1, 1 }, { 1.0f, 0.2f, 1.0f }, 1 },
	{ { -1, -1, -1, 1 }, { 1.0f, 1.0f, 0.2f }, 1 },
	{ {  1, -1, -1, 1 }, { 0.2f, 1.0f, 1.0f }, 1 },
	{ {  1,  1, -1, 1 }, { 1.0f, 0.3f, 0.3f }, 1 },
	{ { -1,  1, -1, 1 }, { 0.2f, 1.0f, 0.3f }, 1 },
};

// 无论画什么图形，都是以三角形逼近的
void draw_plane(device_t* device, int a, int b, int c, int d) {
	vertex_t pos1 = mesh[a], pos2 = mesh[b], pos3 = mesh[c], pos4 = mesh[d];
	device_draw_primitive(device, &pos1, &pos2, &pos3);
	device_draw_primitive(device, &pos3, &pos4, &pos1);

}

// 以 draw_plane 为基础绘制六个面围成体
void draw_box(device_t* device, float theta) {
	matrix_t m;
	matrix_set_rotate(&m, -1, -0.5, 1, theta);
	device->transform.world = m;
	transform_update(&device->transform);

	draw_plane(device, 0, 1, 2, 3);
	draw_plane(device, 7, 6, 5, 4);
	draw_plane(device, 0, 4, 5, 1);
	draw_plane(device, 1, 5, 6, 2);
	draw_plane(device, 2, 6, 7, 3);
	draw_plane(device, 3, 7, 4, 0);
}

void camera_at_zero(device_t* device, float x, float y, float z) {
	point_t eye = { x, y, z, 1 }, at = { 0, 0, 0, 1 }, up = { 0, 0, 1, 1 };
	matrix_set_lookat(&device->transform.view, &eye, &at, &up);
	transform_update(&device->transform);
}

void draw_sphere(device_t* device, vector<vertex_t> vector, float theta) {
	matrix_t m;
	matrix_set_rotate(&m, -1, -0.5, 1, theta);
	device->transform.world = m;
	transform_update(&device->transform);
	int y_segment = 50;
	int x_segment = 50;
	for (int i = 0; i < y_segment; i++)
	{
		for (int j = 0; j < x_segment; j++)
		{
			
			device_draw_primitive(device, &vector[i * (x_segment + 1) + j], &vector[(i+1) * (x_segment + 1) + j], &vector[(i + 1) * (x_segment + 1) + j + 1]);
			device_draw_primitive(device, &vector[i * (x_segment + 1) + j], &vector[(i + 1) * (x_segment + 1) + j + 1], &vector[i * (x_segment + 1) + j + 1]);
		}
	}
}

int main(void)
{
	// 初始化球顶点
	vector<vertex_t> sphere_vertex;
	int x_segment = 50;
	int y_segment = 50;
	for (int y = 0; y <= y_segment; y++) {
		for (int x = 0; x <= x_segment; x++) {
			float xSegment = (float)x / x_segment;
			float ySegment = (float)y / y_segment;
			float x_pos = cos(xSegment * 2.0f * 3.1415926f) * sin(ySegment * 3.1415926f);
			float y_Pos = cos(ySegment * 3.1415926f);
			float z_Pos = sin(xSegment * 2.0f * 3.1415926f) * sin(ySegment * 3.1415926f);			sphere_vertex.push_back({				{x_pos, y_Pos, z_Pos, 1.0f},				{ 1.0f, 0.2f, 0.2f },				1				});		}
	}
	device_t device;
	int states[] = { RENDER_STATE_COLOR, RENDER_STATE_WIREFRAME };
	int indicator = 0;
	int kbhit = 0;
	float alpha = 1;
	float pos = 3.5;
	if (screen_init(800, 600, "3DRender"))
		return -1;
	device_init(&device, 800, 600, screen_fb);
	device.render_state = RENDER_STATE_WIREFRAME;

	camera_at_zero(&device, 3, 0, 0);
	cout << "渲染立方体请输入1，渲染球体请输入2：" << endl;
	int choice;
	cin >> choice;

	while (screen_exit == 0 && screen_keys[VK_ESCAPE] == 0) {
		screen_dispatch();
		device_clear(&device, 1);
		camera_at_zero(&device, pos, 0, 0);

		if (screen_keys[VK_UP]) pos -= 0.01f;
		if (screen_keys[VK_DOWN]) pos += 0.01f;
		if (screen_keys[VK_LEFT]) alpha += 0.01f;
		if (screen_keys[VK_RIGHT]) alpha -= 0.01f;

		if (screen_keys[VK_SPACE]) {
			if (kbhit == 0) {
				kbhit = 1;
				if (++indicator >= 2) indicator = 0;
				device.render_state = states[indicator];
			}
		}
		else {
			kbhit = 0;
		}

		if (choice == 1) {
			draw_box(&device, alpha);
		}
		else if (choice == 2) {
			draw_sphere(&device, sphere_vertex, alpha);
		}
		else {
			return -1;
		}


		screen_update();
		Sleep(1);
	}
	return 0;
}

