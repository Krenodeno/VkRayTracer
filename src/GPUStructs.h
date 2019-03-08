#ifndef GPU_STRUCTS_H
#define GPU_STRUCTS_H

#include <iostream>

#include "bvh.h"
#include "ray.h"
#include "vec.h"

struct RayGPU {
	alignas(16) vec3 origin;
	alignas(16) vec3 direction;
	float tmax;

	RayGPU() : origin(), direction(), tmax(-1.f) {}
	RayGPU(const Point& o, const Point& e) : origin(o), direction(Vector(o, e)), tmax(1.f) {}
	RayGPU(const Point& o, const Vector& d) : origin(o), direction(d), tmax(std::numeric_limits<float>::max()) {}
};

struct HitGPU {
	int id;
	float t, u, v;

	HitGPU() : id(-1), t(0.f), u(0.f), v(0.f) {}

	HitGPU(int id, float t, float u, float v) : id(id), t(t), u(u), v(v) {}

	operator bool() { return id != -1; }
};

struct TriangleGPU {
	alignas(16) vec3 pos;
	alignas(16) vec3 e1;
	alignas(16) vec3 e2;
	int id;

	TriangleGPU(Triangle tri) : pos(tri.p), e1(tri.e1), e2(tri.e2), id(tri.id) {}

	/* Ray - Triangle Intersection */
	HitGPU intersect( const RayGPU &ray, const float htmax ) const
	{
		Vector pvec = cross(Vector(ray.direction), Vector(e2));
		float det = dot(e1, pvec);

		float inv_det = 1 / det;
		Vector tvec(pos, ray.origin);

		float u = dot(tvec, pvec) * inv_det;
		if (u < 0 || u > 1) return HitGPU();

		Vector qvec = cross(tvec, e1);
		float v = dot(ray.direction, qvec) * inv_det;
		if (v < 0 || u + v > 1) return HitGPU();

		float t = dot(e2, qvec) * inv_det;
		if (t > htmax || t < 0) return HitGPU();

		return HitGPU(id, t, u, v);	// p(u, v)= (1 - u - v) * a + u * b + v * c
	}
};

struct BNodeGPU {
	vec3 pmin;
	int next;
	vec3 pmax;
	int skip;

	bool leaf() const { return next < 0; }

	bool node() const { return !leaf(); }
	/** First triangle of the leaf */
	int begin() const { assert(leaf()); return ( -next >> 4 ); }
	/** Last triangle of the leaf */
	int end() const { assert(leaf()); return ( (-next >> 4) + (-next & 0xF) ); }

	/* Ray - Box Intersection */
	NodeHit intersect( const RayGPU& ray, const Vector& invd, const float htmax ) const
	{
		vec3 rmin = pmin;
		vec3 rmax = pmax;
		if(ray.direction.x < 0) std::swap(rmin.x, rmax.x);
		if(ray.direction.y < 0) std::swap(rmin.y, rmax.y);
		if(ray.direction.z < 0) std::swap(rmin.z, rmax.z);
		Vector dmin = (Point(rmin) - Point(ray.origin)) * invd;
		Vector dmax = (Point(rmax) - Point(ray.origin)) * invd;

		float tmin = std::max(dmin.z, std::max(dmin.y, std::max(dmin.x, 0.f)));
		float tmax = std::min(dmax.z, std::min(dmax.y, std::min(dmax.x, htmax)));
		return NodeHit(tmin, tmax);
	}

	BNodeGPU() : skip(-1) {}

};


std::vector<BNodeGPU> transform(const BVH& bvh) {

	assert(!bvh.nodes[bvh.root].leaf());

	std::vector<BNodeGPU> remap;
	remap.reserve(bvh.nodes.size());

	std::vector<int> stackG;
	int current = 0;	// current node in the remapped bvh

	int id = bvh.root;	// current node in the original bvh
	bool cont = true;

	std::vector<int> stack;
	stack.reserve(bvh.nodes.size());

	// First step : root
	{
		const Node& node = bvh.nodes[bvh.root];
		BNodeGPU root;
		root.pmin = vec3(node.bounds.pmin.x, node.bounds.pmin.y, node.bounds.pmin.z);
		root.pmax = vec3(node.bounds.pmax.x, node.bounds.pmax.y, node.bounds.pmax.z);
		root.skip = -1;
		root.next = current + 1;	// doesn't work is root is a leaf

		id = node.left;
		stackG.push_back(current++);
		stack.push_back(bvh.root);
		remap.push_back(std::move(root));
	}

	while (cont) {
		const Node& node = bvh.nodes[id];
		{
			BNodeGPU res;
			res.pmin = vec3(node.bounds.pmin.x, node.bounds.pmin.y, node.bounds.pmin.z);
			res.pmax = vec3(node.bounds.pmax.x, node.bounds.pmax.y, node.bounds.pmax.z);

			if (node.leaf()) {
				assert(node.n() < 15);
				res.next = -( (node.begin() << 4) | node.n() );	// Triangle indices (the offset must be less than 16)
				assert(res.leaf());
			}
			else {
				res.next = current + 1;
				assert(res.node());
			}
			remap.push_back(res);
		}

		// Tree walking :
		// Leaf
		if (node.leaf()) {

			remap[current].skip = remap.size();
			// Right leaf : write skip for each of its parents
			while (bvh.nodes[stack.back()].right == id) {
				//go up
				id = stack.back();
				stack.pop_back();
				current = stackG.back();
				stackG.pop_back();
				//set skip
				remap[current].skip = remap.size();
				if (stack.empty()) break;	// End
			}
			if (stack.empty()) break;
			id = bvh.nodes[stack.back()].right;
			current = remap.size();

		}
		// Not a leaf
		else {
			stack.push_back(id);
			id = node.left;
			stackG.push_back(current++);
		}

		// end of the loop
		if (stack.empty())
			cont = false;
	}

	// little cleanup : set node.skip to null if it is set to bvh.size()
	int cpt = 0;
	for (auto& node : remap) {
		if (node.skip >= remap.size())
			node.skip = -1;
	}

	return remap;
}


HitGPU intersect(const std::vector<BNodeGPU>& pbvh, const std::vector<TriangleGPU>& tris, RayGPU ray) {
	assert(!pbvh.empty());
	assert(!tris.empty());

	Vector invd = Vector(1.f / ray.direction.x, 1.f / ray.direction.y, 1.f / ray.direction.z);

	int id = 0;
	HitGPU hit;
	hit.t = ray.tmax;

	while(id != -1) {
		assert(size_t(id) < pbvh.size());
		const BNodeGPU& node = pbvh[id];
		if (node.leaf()) {
			for (int i = node.begin(); i != node.end(); i++) {
				// test the ray, return the nearest
				HitGPU h = tris[i].intersect(ray, hit.t);
				if (h) {
					hit = h;
				}
			}
			id = node.skip;
		}
		else {
			if (node.intersect(ray, invd, hit.t)) {
				id = node.next;
			}
			else {
				id = node.skip;
			}
		}
	}
	return hit;
}

#endif
