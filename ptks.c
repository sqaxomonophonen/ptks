
#if __STDC_VERSION__ >= 199901L
#define _XOPEN_SOURCE 600
#else
#define _XOPEN_SOURCE 500
#endif /* __STDC_VERSION__ */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>

#include <SDL.h>

#ifndef M_PI
#define M_PI (3.141592653589793)
#endif

struct vec3 {
	float s[3];
};

struct basis {
	struct vec3 u, v, w;
};

// WORLD
struct wall {
	float x0, y0, x1, y1;
	int m0;
	//int portal_wall;
	int portal_sector;

	// derived
	float length;
	struct basis basis;
};

struct sector {
	int wall0;
	int nwalls;
	float z0, z1;
	int m0, m1;

	// derived
	struct basis basis0, basis1;
};


struct worldlet {
	struct wall* walls;
	struct sector* sectors;
	int nsectors;
};


// your random number god
struct rng {
	uint32_t z;
	uint32_t w;
};

static inline uint32_t rng_uint32(struct rng* rng)
{
	/*
	   The MWC generator concatenates two 16-bit multiply-
	   with-carry generators, x(n)=36969x(n-1)+carry,
	   y(n)=18000y(n-1)+carry mod 2^16, has period about
	   2^60 and seems to pass all tests of randomness. A
	   favorite stand-alone generator---faster than KISS,
	   which contains it.
	*/
	rng->z = 36969 * (rng->z & 65535) + (rng->z>>16);
	rng->w = 18000 * (rng->w & 65535) + (rng->w>>16);
	return (rng->z<<16) + rng->w;
}

static inline float rng_float(struct rng* rng)
{
	union {
		uint32_t i;
		float f;
	} magick;
	uint32_t r = rng_uint32(rng);
	magick.i = (r & 0x007fffff) | (127 << 23);
	return magick.f - 1;
}

static inline void rng_seed(struct rng* rng, uint32_t seed)
{
	rng->z = 654654 + seed;
	rng->w = 7653234 + seed * 69069;
}


////

static void sdl_panic()
{
	fprintf(stderr, "SDL: %s\n", SDL_GetError());
	exit(EXIT_FAILURE);
}

static inline void vec3_sub(struct vec3* dst, struct vec3* a, struct vec3* b)
{
	for (int i = 0; i < 3; i++) {
		dst->s[i] = a->s[i] - b->s[i];
	}
}

static inline void vec3_cross(struct vec3* restrict dst, struct vec3* restrict u, struct vec3* restrict v)
{
	float s[3];
	for (int i = 0; i < 3; i++) {
		int j = (i + 1) % 3;
		int k = (i + 2) % 3;
		s[i] = u->s[j] * v->s[k] - u->s[k] * v->s[j];
	}
	for (int i = 0; i < 3; i++) {
		dst->s[i] = s[i];
	}
}

static inline void vec3_add_scaled_inplace(struct vec3* restrict dst, struct vec3* restrict src, float s)
{
	for (int i = 0; i < 3; i++) {
		dst->s[i] += src->s[i] * s;
	}
}

static inline void vec3_scale(struct vec3* restrict v, float s)
{
	for (int i = 0; i < 3; i++) {
		v->s[i] *= s;
	}
}

static inline void vec3_normalize(struct vec3* restrict v)
{
	float d2 = 0;
	for (int i = 0; i < 3; i++) {
		d2 += v->s[i] * v->s[i];
	}
	float rd = 1 / sqrtf(d2);
	for (int i = 0; i < 3; i++) {
		v->s[i] *= rd;
	}
}

static inline void vec3_mult_inplace(struct vec3* restrict u, struct vec3* restrict v)
{
	for (int i = 0; i < 3; i++) {
		u->s[i] *= v->s[i];
	}
}

static inline void vec3_add_inplace(struct vec3* restrict u, struct vec3* restrict v)
{
	for (int i = 0; i < 3; i++) {
		u->s[i] += v->s[i];
	}
}

static inline int vec3_zero(struct vec3* restrict v)
{
	return v->s[0] == 0 && v->s[1] == 0 && v->s[2] == 0;
}

static inline float vec3_dot(struct vec3* a, struct vec3* b)
{
	return a->s[0] * b->s[0] + a->s[1] * b->s[1];
}

static inline float vec3_length(struct vec3* v)
{
	return sqrtf(vec3_dot(v, v));
}

struct camera {
	struct vec3 origin;
	struct vec3 front;
	struct vec3 right;
	struct vec3 up;
};

static void camera_point_at(struct camera* c, struct vec3* p, struct vec3* up)
{
	vec3_sub(&c->front, p, &c->origin);
	vec3_normalize(&c->front);
	vec3_cross(&c->right, &c->front, up);
	vec3_cross(&c->up, &c->front, &c->right);
}


struct worldlet* worldlet_0()
{
	/*
	returns this!
	    0---1
	    |   |
	    |   |
	    |   |
	2---3---4---5
	|.         .|
	| .       . |
	|  .     .  |
	6---6---7---7

	0: (1,0)
	1: (2,0)
	2: (0,1)
	3: (1,1)
	4: (2,1)
	5: (3,1)
	6: (0,2)
	7: (2,2)

	*/
	struct worldlet* worldlet = malloc(sizeof(struct worldlet));
	int nwalls = 4 + 6;
	worldlet->walls = malloc(nwalls * sizeof(struct wall));

	int points[8][2] = {
		{1,0},
		{2,0},
		{0,1},
		{1,1},
		{2,1},
		{3,1},
		{0,2}, //{1,2},
		{3,2} //{2,2}
	};

	int w0[5] = {0,1,4,3,0};
	int m00[4] = {2,1,-1,1};
	int w1[7] = {2,3,4,5,7,6,2};
	int m01[6] = {1,-1,1,1,1,1};

	for (int i = 0; i < 4; i++) {
		struct wall* w = &worldlet->walls[i];
		w->x0 = points[w0[i]][0];
		w->y0 = points[w0[i]][1];
		w->x1 = points[w0[i+1]][0];
		w->y1 = points[w0[i+1]][1];
		w->m0 = m00[i];
		//w->portal_wall = -1;
		w->portal_sector = -1;
	}

	for (int i = 0; i < 6; i++) {
		struct wall* w = &worldlet->walls[i+4];
		w->x0 = points[w1[i]][0];
		w->y0 = points[w1[i]][1];
		w->x1 = points[w1[i+1]][0];
		w->y1 = points[w1[i+1]][1];
		w->m0 = m01[i];
		//w->portal_wall = -1;
		w->portal_sector = -1;
	}

	worldlet->nsectors = 2;
	worldlet->sectors = malloc(worldlet->nsectors * sizeof(struct sector));

	float ceil = 1;
	float floor = 0;

	{
		struct sector* s = &worldlet->sectors[0];
		s->wall0 = 0;
		s->nwalls = 4;
		s->m0 = 0;
		s->m1 = 0;
		s->z0 = floor;
		s->z1 = ceil;
	}
	{
		struct sector* s = &worldlet->sectors[1];
		s->wall0 = 4;
		s->nwalls = 6;
		s->m0 = 0;
		s->m1 = 0;
		s->z0 = floor;
		s->z1 = ceil;
	}
	{
		struct wall* w = &worldlet->walls[2];
		w->portal_sector = 1;
		//w->portal_wall = 5;
	}
	{
		struct wall* w = &worldlet->walls[5];
		w->portal_sector = 0;
		//w->portal_wall = 2;
	}

	return worldlet;
}


static inline void basis_for_normal(struct basis* restrict basis, struct vec3* restrict normal)
{
	basis->w = *normal;
	vec3_normalize(&basis->w);

	struct vec3 bx = {{1,0,0}};
	struct vec3 by = {{0,1,0}};
	vec3_cross(&basis->u, &basis->w, fabs(basis->w.s[0]) > 0.2 ? &by : &bx);
	vec3_normalize(&basis->u);

	vec3_cross(&basis->v, &basis->w, &basis->u);
	vec3_normalize(&basis->v);
}

static void worldlet_update(struct worldlet* w)
{
	for (int i = 0; i < w->nsectors; i++) {
		struct sector* sector = &w->sectors[i];

		struct vec3 normal0 = {{0,0,1}};
		struct vec3 normal1 = {{0,0,-1}};
		basis_for_normal(&sector->basis0, &normal0);
		basis_for_normal(&sector->basis1, &normal1);

		for (int j = sector->wall0; j < (sector->wall0 + sector->nwalls); j++) {
			struct wall* wall = &w->walls[j];
			float dx = wall->x1 - wall->x0;
			float dy = wall->y1 - wall->y0;
			wall->length = sqrtf(dx*dx + dy*dy);

			struct vec3 normal = {{-dy, dx, 0}};
			basis_for_normal(&wall->basis, &normal);
		}
	}
}

#define COS_DIST_SIZE (1<<10)

struct renderer {
	uint8_t* materials;

	struct worldlet* worldlet;
	struct camera* camera;
	int sector0;
	int width;
	int height;

	struct vec3 dbase;
	struct vec3 dx;
	struct vec3 dy;

	struct vec3* cos_dist;
};

static void renderer_init(struct renderer* r)
{

	int fd = open("mat.dat", O_RDONLY);
	if (fd == -1) {
		perror("mat.dat");
		exit(1);
	}

	struct stat st;
	if (fstat(fd, &st) == -1) {
		perror("mat.dat");
		exit(1);
	}

	r->materials = mmap(NULL, st.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
	if (r->materials == MAP_FAILED) {
		perror("mmap");
		exit(1);
	}

	r->cos_dist = malloc(COS_DIST_SIZE * sizeof(struct vec3));
	struct rng rng;
	rng_seed(&rng, 6667);
	for (int i = 0; i < COS_DIST_SIZE; i++) {
		float r1 = rng_float(&rng) * 2 * M_PI;
		float r2 = rng_float(&rng);
		float r2s = sqrtf(r2);

		struct vec3 c = {{
			cosf(r1) * r2s,
			sinf(r1) * r2s,
			sqrtf(1-r2)
		}};

		r->cos_dist[i] = c;
	}
}

static void renderer_reset(struct renderer* r, struct worldlet* worldlet, struct camera* camera, int sector0, int width, int height)
{
	r->worldlet = worldlet;
	r->camera = camera;
	r->sector0 = sector0;
	r->width = width;
	r->height = height;

	r->dbase = camera->front;
	float viewplane_distance = 0.3f;
	vec3_scale(&r->dbase, viewplane_distance);

	vec3_add_scaled_inplace(&r->dbase, &camera->right, -0.5f);
	vec3_add_scaled_inplace(&r->dbase, &camera->up, (float)r->height/(float)r->width/2.0f);

	r->dx = camera->right;
	vec3_scale(&r->dx, 1.0f / (float)r->width);

	r->dy = camera->up;
	vec3_scale(&r->dy, -1.0f / (float)r->width);
}

struct hit {
	struct vec3 p0;
	int material;
	float u, v;
	struct basis* basis;
};

inline static void renderer_sector_trace(
	struct renderer* restrict r,
	struct vec3* restrict origin,
	struct vec3* restrict dir,
	struct sector** restrict sectorp,
	struct hit* restrict hit)
{
	// shut up various compiler warnings:
	struct vec3 zero = {{0,0,0}};
	hit->p0 = zero;

	for (;;) {

		struct sector* sector = *sectorp;

		float min_d = INFINITY;

		int next_sector = -1;

		float dz = dir->s[2];
		float l0 = origin->s[2]; // point on line
		if (dz < 0) { // down ray; trace hit against floor
			float p0 = sector->z0; // point on floor
			float d = (p0 - l0) / dz;
			if (d < min_d) {
				min_d = d;
				struct vec3 p0 = *origin;
				vec3_add_scaled_inplace(&p0, dir, d);
				hit->p0 = p0;
				hit->basis = &sector->basis0;
				hit->material = sector->m0;
				hit->u = p0.s[0];
				hit->v = p0.s[1];
			}
		} else { // up ray; trace hit against ceiling
			float p0 = sector->z1; // point on ceiling
			float d = (p0 - l0) / dz;
			if (d < min_d) {
				min_d = d;
				struct vec3 p0 = *origin;
				vec3_add_scaled_inplace(&p0, dir, d);
				hit->p0 = p0;
				hit->basis = &sector->basis1;
				hit->material = sector->m1;
				hit->u = p0.s[0];
				hit->v = p0.s[1];
			}
		}

		// trace hits against walls
		for (int i = 0; i < sector->nwalls; i++) {
			int wall_index = sector->wall0 + i;

			struct wall* wall = &r->worldlet->walls[wall_index];

			float ddx = dir->s[0];
			float ddy = dir->s[1];

			float wdx = wall->x1 - wall->x0;
			float wdy = wall->y1 - wall->y0;

			float wnx = -wdy;
			float wny = wdx;

			float dot = ddx*wnx + ddy*wny;

			if (dot >= 0.0f) continue;

			/*
			ray : p + t*r
			wall : q + u*s

			t = ((q-p) x s) / (r x s)
			u = ((q-p) x r) / (r x s)

			(cross product v x w = vx * wy - vy * wx)
			*/

			// rxs
			float rxs1 = 1.0f / (ddx * wdy - ddy * wdx);

			// q-p
			float dx = wall->x0 - origin->s[0];
			float dy = wall->y0 - origin->s[1];

			// ((q-p) x s) / (r x s)
			float t = (dx*wdy - dy*wdx) * rxs1;

			if (t >= min_d) continue;

			// ((q-p) x r) / (r x s)
			float u = (dx*ddy - dy*ddx) * rxs1;

			if (u < 0.0f || u > 1.0f) continue; // outside line segment; disregard

			min_d = t;
			struct vec3 p0 = *origin;
			vec3_add_scaled_inplace(&p0, dir, t);
			hit->p0 = p0;
			hit->basis = &wall->basis;
			hit->material = wall->m0;
			hit->u = u * wall->length;
			hit->v = p0.s[2];

			next_sector = wall->portal_sector;
		}

		if (next_sector == -1) break;

		(*sectorp) = &r->worldlet->sectors[next_sector];
		hit->material = 0; // fix for rare case where a ray goes through a portal and hits nothing
	}
}

#define TEXSIZE (128)

inline static void material_lookup(struct renderer* restrict r, struct hit* restrict hit, struct vec3* restrict color, struct vec3* restrict emission)
{
	int u = (int)(hit->u * (float)TEXSIZE) & (TEXSIZE-1);
	int v = (int)(hit->v * (float)TEXSIZE) & (TEXSIZE-1);
	int m = hit->material;

	float scale = 1.0f / 256.0f;

	uint8_t* texel = r->materials + m * (TEXSIZE*TEXSIZE*6) + u * 6 + v * TEXSIZE * 6;
	for (int i = 0; i < 3; i++) {
		color->s[i] = (float)texel[i] * scale;
		emission->s[i] = (float)texel[i + 3] * scale;
	}
}

inline static void renderer_sample_pixel(struct renderer* restrict r, struct vec3* restrict pixel, float x, float y, struct rng* restrict rng)
{
	struct vec3 o = r->camera->origin;
	struct vec3 d = r->dbase;
	vec3_add_scaled_inplace(&d, &r->dx, x);
	vec3_add_scaled_inplace(&d, &r->dy, y);

	struct sector* sector = &r->worldlet->sectors[r->sector0];
	struct hit hit;
	hit.material = 0;
	hit.u = 0;
	hit.v = 0;
	struct basis fake_basis;
	hit.basis = &fake_basis;
	int depth = 0;

	struct vec3 diffuse = {{1,1,1}};

	for (;;) {

		renderer_sector_trace(r, &o, &d, &sector, &hit);
		o = hit.p0;

		struct vec3 color;
		struct vec3 emission;

		material_lookup(r, &hit, &color, &emission);

		if (!vec3_zero(&emission)) {
			struct vec3 contrib;
			contrib = emission;
			vec3_mult_inplace(&contrib, &diffuse);
			vec3_add_inplace(pixel, &contrib);
		}

		if (vec3_zero(&color)) {
			// no diffuse; we're done
			return;
		}

		vec3_mult_inplace(&diffuse, &color);

		// russian roulette
		if (depth > 3) {
			return; // bang; you're dead! (XXX FIXME)
			float p =
				color.s[0] > color.s[1] && color.s[0] > color.s[2] ? color.s[0] :
				color.s[1] > color.s[2] ? color.s[1] :
				color.s[2];
			if (rng_float(rng) < p) {
				vec3_scale(&color, 1/p);
			} else {
				// kill "recursion"
				return;
			}
		}

		// find new cosine-distributed direction. or something.
		struct vec3* rr = &r->cos_dist[rng_uint32(rng) & (COS_DIST_SIZE-1)];
		struct vec3 new_d = {{0,0,0}};
		vec3_add_scaled_inplace(&new_d, &hit.basis->u, rr->s[0]);
		vec3_add_scaled_inplace(&new_d, &hit.basis->v, rr->s[1]);
		vec3_add_scaled_inplace(&new_d, &hit.basis->w, rr->s[2]);

		d = new_d;

		depth++;

	}
}


static inline int clamp(int v, int min, int max)
{
	return v > max ? max : v < min ? min : v;
}

static inline uint8_t f2c(float f)
{
	return clamp(f * 255, 0, 255);
}





int worldlet_find_sector(struct worldlet* w, struct vec3* o)
{
	for (int i = 0; i < w->nsectors; i++) {
		struct sector* sector = &w->sectors[i];

		int intersect_count = 0;

		for (int j = 0; j < sector->nwalls; j++) {
			int wall_index = sector->wall0 + j;

			struct wall* wall = &w->walls[wall_index];

			float wdx = wall->x1 - wall->x0;
			float wdy = wall->y1 - wall->y0;

			if (wdy == 0) continue;

			// rxs
			float rxs1 = 1.0f / wdy;

			// q-p
			float dx = wall->x0 - o->s[0];
			float dy = wall->y0 - o->s[1];

			// ((q-p) x s) / (r x s)
			float t = (dx*wdy - dy*wdx) * rxs1;

			// ((q-p) x r) / (r x s)
			float u = -dy * rxs1;

			if (t < 0 || u <= 0.0f || u > 1.0f) continue;

			intersect_count++;
		}

		// if ray intersects with an odd number of ecges, we're inside the sector
		if (intersect_count & 1) return i;
	}
	return -1;
}

static uint8_t blend(uint8_t a, uint8_t b)
{
	return (a + b) >> 1;
}

int main(int argc, char** argv)
{
	if(SDL_Init(SDL_INIT_VIDEO) != 0) sdl_panic();

	struct timespec t0;
	struct timespec t1;

	int anti_alias_root = 7;
	int anti_alias_squared = anti_alias_root * anti_alias_root;
	int width = 192;
	int height = 108;
	int nsamples = 1;

	float brightness_scale = 1.0f / ((float)anti_alias_squared) / (float)nsamples;
	float complexity = (float)anti_alias_squared * (float)nsamples * ((float)width * (float)height) / (320.0 * 240.0);

	uint8_t* pixels = malloc(width * height * 4);
	uint8_t* pixels2 = malloc(width * height * 4);
	bzero(pixels2, width * height * 4);

	float iris = 1;

	SDL_Window* sdl_window = SDL_CreateWindow(
		"PTKS",
		SDL_WINDOWPOS_UNDEFINED,
		SDL_WINDOWPOS_UNDEFINED,
		width,
		height,
		0);

	SDL_Renderer* sdl_renderer = SDL_CreateRenderer(sdl_window, -1, 0);

	SDL_Texture* sdl_texture = SDL_CreateTexture(
		sdl_renderer,
		SDL_PIXELFORMAT_ARGB8888,
		SDL_TEXTUREACCESS_STREAMING,
		width,
		height);

	struct worldlet* worldlet = worldlet_0();
	worldlet_update(worldlet);

	int exiting = 0;
	int frame = 0;
	int animate = 1;
	float angle = 0;

	struct renderer renderer;
	renderer_init(&renderer);

	int move_forward = 0;
	int move_backward = 0;
	int turn_left = 0;
	int turn_right = 0;
	struct vec3 co = {{1.5, 1.5, 0.5}}; // ORIGIN
	int sector0 = worldlet_find_sector(worldlet, &co);

	if (sector0 == -1) {
		fprintf(stderr, "initial origin outside world!\n");
		exit(1);
	}

	while (!exiting) {
		clock_gettime(CLOCK_MONOTONIC,&t0);

		struct camera camera;

		float turn_speed = 0.006f * complexity;
		if (turn_right) {
			angle -= turn_speed;
		}
		if (turn_left) {
			angle += turn_speed;
		}
		struct vec3 cdir = {{+sin(angle), -cos(angle), 0}};
		struct vec3 cp = co;
		vec3_add_inplace(&cp, &cdir);
		struct vec3 up = {{0,0,1}};
		float move_speed = 0.003f * complexity;
		if (move_forward) {
			vec3_add_scaled_inplace (&co, &cdir, move_speed);
		}
		if (move_backward) {
			vec3_add_scaled_inplace (&co, &cdir, -move_speed);
		}
		camera.origin = co;
		camera_point_at(&camera, &cp, &up);

		int new_sector = worldlet_find_sector(worldlet, &co);
		if (new_sector != -1) {
			sector0 = new_sector;
		}

		renderer_reset(&renderer, worldlet, &camera, sector0, width, height);
		float anti_alias_step = 1/(float)anti_alias_root;

		double emission = 0;

		#pragma omp parallel for schedule(dynamic, 8)
		for (int y = 0; y < height; y++) {
			struct rng rng;
			rng_seed(&rng, y*69069 + frame);

			int i = y * width;
			for (int x = 0; x < width; x++) {
				struct vec3 pixel = {{0,0,0}};
				for (float ay = 0; ay < 1; ay+=anti_alias_step) {
					for (float ax = 0; ax < 1; ax+=anti_alias_step) {
						for (int s = 0; s < nsamples; s++) {
							renderer_sample_pixel(&renderer, &pixel, x + ax , y + ay, &rng);
						}
					}
				}
				//emission += vec3_length(&pixel);
				emission += pixel.s[0] * 0.2126 + pixel.s[1] * 0.7152 + pixel.s[2] * 0.0722;
				vec3_scale(&pixel, iris);
				uint8_t red = f2c(pixel.s[0] * brightness_scale);
				uint8_t green = f2c(pixel.s[1] * brightness_scale);
				uint8_t blue = f2c(pixel.s[2] * brightness_scale);
				pixels[i * 4] = blend(pixels2[i * 4], blue);
				pixels[i * 4 + 1] = blend(pixels2[i * 4 + 1], green);
				pixels[i * 4 + 2] = blend(pixels2[i * 4 + 2], red);
				pixels[i * 4 + 3] = 255;
				i++;
			}
		}

		double target_factor = 70;
		double target_nfactor = width * height; // * anti_alias_squared * nsamples;
		double target_emission = target_nfactor * target_factor;
		double target_speed = 0.13;

		iris += ((target_emission/emission) - iris) * target_speed;


		memcpy(pixels2, pixels, width * height * 4);

		clock_gettime(CLOCK_MONOTONIC, &t1);

		float dt = (t1.tv_sec - t0.tv_sec) + (float)(t1.tv_nsec - t0.tv_nsec) * 1e-9;
		printf("frame: %f seconds\n", dt);
		printf("emission: %f / corrected: %f\n", emission / target_nfactor, emission * iris / target_nfactor);

		SDL_UpdateTexture(sdl_texture, NULL, pixels, width * sizeof(uint32_t));
		SDL_RenderClear(sdl_renderer);
		SDL_RenderCopy(sdl_renderer, sdl_texture, NULL, NULL);
		SDL_RenderPresent(sdl_renderer);

		if (animate) frame++;

		SDL_Event e;
		while (SDL_PollEvent(&e)) {
			if (e.type == SDL_KEYDOWN) {
				if (e.key.keysym.sym == SDLK_ESCAPE || e.key.keysym.sym == SDLK_q) {
					exiting = 1;
				}
				if (e.key.keysym.sym == SDLK_SPACE) {
					animate = !animate;
				}
				if (e.key.keysym.sym == SDLK_UP) {
					move_forward = 1;
				}
				if (e.key.keysym.sym == SDLK_DOWN) {
					move_backward = 1;
				}
				if (e.key.keysym.sym == SDLK_LEFT) {
					turn_left = 1;
				}
				if (e.key.keysym.sym == SDLK_RIGHT) {
					turn_right = 1;
				}
			}
			if (e.type == SDL_KEYUP) {
				if (e.key.keysym.sym == SDLK_UP) {
					move_forward = 0;
				}
				if (e.key.keysym.sym == SDLK_DOWN) {
					move_backward = 0;
				}
				if (e.key.keysym.sym == SDLK_LEFT) {
					turn_left = 0;
				}
				if (e.key.keysym.sym == SDLK_RIGHT) {
					turn_right = 0;
				}
			}
		}

	}

	return 0;
}



