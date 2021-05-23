#pragma once
#ifndef CCVT_OPTIMIZER_H
#define CCVT_OPTIMIZER_H

#include <algorithm>
#include <assert.h>
#include <limits>
#include <list>
#include <math.h>
#include <map>
#include <list>
#include <queue>
#include <vector>

namespace ccvt {

	template<class Site, class Point, class Metric>
	class Optimizer {

	private:

		struct Entry;

		/* The initialization of the optimizer uses a kd-tree to initially assign
		   the points to the sites. A faster method is to use the ordinary Voronoi
		   tessellation for that matter. However, the kd-tree is used because it
		   can easily be implemented for the n-dimensional case. */
		class KdTree {//KdTree

		public:

			KdTree(std::vector<Entry> &entries, const Metric& metric, const int dimensions)
				: dimensions_(dimensions), root_(NULL), metric_(metric), lastNearestNeighborNode_(NULL) {
				if (!entries.empty()) {
					std::vector<Entry*> entriesPtr(entries.size());
					for (int i = 0; i < static_cast<int>(entriesPtr.size()); ++i) {
						entriesPtr[i] = &entries[i];
					}
					build(root_, entriesPtr, 0, static_cast<int>(entriesPtr.size() - 1), 0);
				}
			}

			~KdTree() {
				if (root_ != NULL) {
					delete root_;
				}
			}

			bool empty() {
				return root_ == NULL || root_->disabled;
			}

			Entry* nearest_neighbor(const Point& point) {
				if (root_ == NULL || root_->disabled) {
					return NULL;
				}
				double minDistance = std::numeric_limits<double>::max();
				NearestNeighborResult result = nearest_neighbor(root_, point, 0, minDistance);
				lastNearestNeighborNode_ = result.first;
				return result.second;
			}

			void disable_last_nearest_neighbor_node() {
				if (lastNearestNeighborNode_ != NULL) {
					lastNearestNeighborNode_->disabled = true;
					Node* parentNode = lastNearestNeighborNode_->parent;
					while (parentNode != NULL && parentNode->left->disabled && parentNode->right->disabled) {
						parentNode->disabled = true;
						parentNode = parentNode->parent;
					}
					lastNearestNeighborNode_ = NULL;
				}
			}

		private:

			KdTree(const KdTree& kdTree);
			KdTree& operator=(const KdTree& kdTree);

			struct Node {

				Node()
					: median(0), left(NULL), right(NULL), entry(NULL), parent(NULL), disabled(false) {
				}

				~Node() {
					if (left != NULL) {
						delete left;
					}
					if (right != NULL) {
						delete right;
					}
				}

				double      median;
				Node*       left;
				Node*       right;
				Entry*      entry;
				Node*       parent;
				bool        disabled;

			};

			typedef std::pair<Node*, Entry*> NearestNeighborResult;

			struct Less {

				Less(int dimension)
					: dimension(dimension) {
				}

				bool operator()(const Entry* p1, const Entry* p2) {
					return p1->site->location[dimension] < p2->site->location[dimension];
				}

				int dimension;

			};

			void build(Node*& node, std::vector<Entry*>& entries, const int min, const int max, const int depth) {
				node = new Node();
				if (min == max) {
					node->entry = entries[min];
					return;
				}
				std::sort(entries.begin() + min, entries.begin() + max + 1, Less(depth % dimensions_));//��Χ��������dimensions 
				int medianIndex = (min + max) / 2;
				node->median = (entries[medianIndex]->site->location[depth % dimensions_] + entries[medianIndex + 1]->site->location[depth % dimensions_]) / 2;
				build(node->left, entries, min, medianIndex, depth + 1);
				build(node->right, entries, medianIndex + 1, max, depth + 1);
				node->left->parent = node;
				node->right->parent = node;
			}

			NearestNeighborResult nearest_neighbor(Node* node, const Point& point, const int depth, double& minDistance) {
				if (node->disabled) {
					return NearestNeighborResult(reinterpret_cast<Node*>(NULL), reinterpret_cast<Entry*>(NULL));
				}

				if (node->entry != NULL) {
					minDistance = metric_.distance(point, node->entry->site->location);
					return NearestNeighborResult(node, node->entry);
				}

				Node* child;
				double otherDistance;
				if (point[depth % dimensions_] <= node->median) {
					child = node->left;
					otherDistance = node->median - point[depth % dimensions_];
				}
				else {
					child = node->right;
					otherDistance = point[depth % dimensions_] - node->median;
				}

				NearestNeighborResult result = nearest_neighbor(child, point, depth + 1, minDistance);

				if (minDistance > otherDistance) {
					double newMinDistance = minDistance;
					Node* otherChild = node->left;
					if (otherChild == child) {
						otherChild = node->right;
					}
					NearestNeighborResult newResult = nearest_neighbor(otherChild, point, depth + 1, newMinDistance);
					if (newMinDistance < minDistance) {
						result = newResult;
						minDistance = newMinDistance;
					}
				}

				return result;
			}

			const int     dimensions_;
			Node*         root_;
			const Metric& metric_;
			Node*         lastNearestNeighborNode_;

		};

	public:

		Optimizer() {
		}

		void clear() {
			entries_.clear();
			mapping_.clear();
			sites_.clear();
			metric_ = Metric();
		}

		void initialize(typename std::list<Site>& sites, typename std::list<Point>& points, const Metric& metric) {
			clear();

			metric_ = metric;

			int sitesSize = static_cast<int>(sites.size());
			sites_.reserve(sitesSize);//锟斤拷锟斤拷站锟斤拷锟斤拷目锟斤拷锟节达拷占锟?锟斤拷锟芥储锟斤拷锟斤拷锟斤拷锟斤拷站锟斤拷 
			entries_.reserve(sitesSize);//
			int sumCapacities = 0;
			for (int i = 0; i < sitesSize; ++i) {
				Site& site = sites.front();
				sumCapacities += site.capacity;
				if (site.capacity > 0) {
					sites_.push_back(site);
					entries_.push_back(Entry(&sites_.back()));
					mapping_.insert(std::make_pair(site.id, &entries_.back()));
				}
				sites.pop_front();
			}
			assert(sumCapacities == points.size()); // the sum of the site capacities must be equal to the number of points
			sites_.resize(sites_.size());//锟斤拷锟斤拷站锟斤拷锟斤拷锟斤拷锟斤拷占锟节达拷锟叫?
			entries_.resize(entries_.size());

			KdTree kdTree(entries_, metric_, Point::D);//锟斤拷锟斤拷Kd锟斤拷 
			while (!points.empty() && !kdTree.empty()) {
				const Point& point = points.back();
				//std::cout << point.x << " " << point.y << std::endl;
				Entry* entry = kdTree.nearest_neighbor(point);
				entry->points.push_back(point);
				if (static_cast<int>(entry->points.size()) == entry->site->capacity) {
					entry->points.resize(entry->points.size());
					kdTree.disable_last_nearest_neighbor_node();
				}
				points.pop_back();
			}

			for (unsigned int i = 0; i < entries_.size(); ++i) {
				entries_[i].energy = 0;
				int pointsSize = static_cast<int>(entries_[i].points.size());
				for (int j = 0; j < pointsSize; ++j) {
					entries_[i].energy += energy(entries_[i].points[j], entries_[i].site);
				}
				entries_[i].update(metric_);
			}
		}

		bool optimize(const bool centroidalize) {
			int entriesSize = static_cast<int>(entries_.size());
			std::vector<bool> stability(entriesSize, true);
			for (int i = 0; i < entriesSize; ++i) {
				for (int j = i + 1; j < entriesSize; ++j) {
					Entry* entry1 = &entries_[i];
					Entry* entry2 = &entries_[j];

//					if (entry1->stable && entry2->stable) {
//						std::cout << "OK-" << std::endl;
//					}

					if (entry1->stable && entry2->stable ||
						metric_.distance(entry1->bounding.center, entry2->bounding.center) > entry1->bounding.radius + entry2->bounding.radius) {
						continue;
					}

					if (entry1->points.size() > entry2->points.size()) {
						std::swap(entry1, entry2);
					}

					Site* site1 = entry1->site;
					Site* site2 = entry2->site;
					//std::cout << site1->location.x << " " << site2->location.x << std::endl;
					std::vector<Point>* points1 = &entry1->points;
					std::vector<Point>* points2 = &entry2->points;

					double maxSquaredRadius = std::max(entry1->bounding.squaredRadius, entry2->bounding.squaredRadius);
					//printf("%d %d %f\n", i, j, maxSquaredRadius);

					typename Candidate::Vector candidates1(points1->size());
					int size = static_cast<int>(points1->size());
					int count = 0;
					double testE = 0;
					for (int k = 0; k < size; ++k) {
						Point& point = (*points1)[k];
						if (metric_.distance_square(point, entry2->bounding.center) <= maxSquaredRadius) {
							candidates1[count++] = Candidate(&point, energy(point, site1), energy(point, site2));
							//forecast the next location and calculate the max deltaE
							Point p(point.x, point.y);
							testE = std::max(testE, energy(p, site1) - energy(p, site2));
						}
					}
					if (count == 0) {
						continue;
					}
					candidates1.resize(count);
					std::make_heap(candidates1.begin(), candidates1.end());
					//std::cout << "OK1"<<std::endl;

					double minEnergy = -(candidates1.front().energySelf - candidates1.front().energyOther);
					//std::cout << minEnergy << std::endl;
					typename Candidate::Vector candidates2(points2->size());
					size = static_cast<int>(points2->size());
					count = 0;
					for (int k = 0; k < size; ++k) {
						Point& point = (*points2)[k];
						if (metric_.distance_square(point, entry1->bounding.center) <= maxSquaredRadius) {
							double eSelf = energy(point, site2);
							double eOther = energy(point, site1);
							if (eSelf - eOther > minEnergy) {
								candidates2[count++] = Candidate(&point, eSelf, eOther);
								Point p(point.x, point.y);
								testE = std::max(testE, eSelf - eOther);
							}
						}
					}
					if (count == 0) {
						continue;
					}
					candidates2.resize(count);
					std::make_heap(candidates2.begin(), candidates2.end());
					//std::cout << "OK2" << std::endl;

					int maxSwaps = static_cast<int>(std::min(candidates1.size(), candidates2.size()));
					int swaps;
					for (swaps = 0; swaps < maxSwaps; ++swaps) {
						Candidate& candidate1 = candidates1.front();
						Candidate& candidate2 = candidates2.front();
						if (candidate1.energySelf - candidate1.energyOther + candidate2.energySelf - candidate2.energyOther <= 0) {
							break;
						}
						std::swap(*candidate1.point, *candidate2.point);
						entry1->energy += candidate2.energyOther - candidate1.energySelf;
						entry2->energy += candidate1.energyOther - candidate2.energySelf;
						std::pop_heap(candidates1.begin(), candidates1.end() - swaps);
						std::pop_heap(candidates2.begin(), candidates2.end() - swaps);
					}

					if (swaps > 0) {
						num_swap += swaps;
						stability[i] = false;
						stability[j] = false;
						if (centroidalize) {
							entry1->site->location = metric_.centroid(entry1->site->location, entry1->points);
							entry2->site->location = metric_.centroid(entry2->site->location, entry2->points);
						}
						entry1->update(metric_);
						entry2->update(metric_);
					}
				}
			}

			bool stable = true;
			for (int i = 0; i < entriesSize; ++i) {
				entries_[i].stable = stability[i];
				stable &= stability[i];
			}
			return stable;
		}

		double energy() const {
			double e = 0;
			int entriesSize = static_cast<int>(entries_.size());
			for (int i = 0; i < entriesSize; ++i) {
				e += entries_[i].energy;
			}
			return e;
		}

		inline double site_energy(const int id) const {
			typename Entry::MapPtr::const_iterator it = mapping_.find(id);
			if (it == mapping_.end()) {
				return 0;
			}
			return it->second->energy;
		}

		inline const std::vector<Point>* site_points(const int id) const {
			typename Entry::MapPtr::const_iterator it = mapping_.find(id);
			if (it == mapping_.end()) {
				return NULL;
			}
			return &it->second->points;
		}

		inline bool site_stable(const int id) const {
			typename Entry::MapPtr::const_iterator it = mapping_.find(id);
			return it == mapping_.end() || it->second->stable;
		}

		inline const std::vector<Site>& sites() {
			return sites_;
		}

		inline const double x_of_sites(int id) {
			return entries_[id].site->location.x;
		}

		inline const double y_of_sites(int id) {
			return entries_[id].site->location.y;
		}

		inline const std::vector<Point>& entryPoint(int id) {
			return entries_[id].points;
		}

		inline void update_site_location(const int id, const Point& location) {
			typename Entry::MapPtr::const_iterator it = mapping_.find(id);
			typename Entry::MapPtr::const_iterator it = mapping_.find(id);
			if (it != mapping_.end()) {
				it->second->site->location = location;
				it->second->update(metric_);
			}
		}

		inline void updatePointFather() {
			int entriesSize = static_cast<int>(entries_.size());
			for (int i = 0; i < entriesSize; i++) {
				Entry* entry = &entries_[i];
				std::vector<Point>* points = &entry->points;
				int size = static_cast<int>(points->size());
				for (int j = 0; j < size; j++) {
					Point& point = (*points)[j];
					point.changeFather(entry->site -> id);
					//std::cout << (*points)[j].getFather() << std::endl;
				}
			}
		}

		inline void updatePointLocation() {
			int entriesSize = static_cast<int>(entries_.size());
			for (int i = 0; i < entriesSize; i++) {
				Entry* entry = &entries_[i];
				entry->stable = false;
				std::vector<Point>* points = &entry->points;
				int size = static_cast<int>(points->size());
				for (int j = 0; j < size; j++) {
					Point& point = (*points)[j];
					//std::cout << point.x << " " << point.y << "||";
					point.update();
					checkBox(point,10.0/100);
					//std::cout << point.father << std::endl;
					//std::cout << point.x << " " << point.y << std::endl;
				}
			}
		}

		inline int getSwapNumber() {
			int num = num_swap;
			num_swap = 0;
			return num;
		}

		inline void initializeGird(const int n,const double L) {
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					Grids.push_back(Grid(i, j, L));
				}
			}
		}

		inline void clear2() {
			for (int i = 0; i < Grids.size(); i++) {
				std::vector<std::vector<Point2>>().swap(Grids[i].points);
			}
			std::priority_queue<SitePair>().swap(sitePairs);
			std::vector<SitePairPoints>().swap(sitePairsPoints);
		}

		inline void checkBox(const Point2 p,const double L) {
			int n = Grids.size();
			n = sqrt(n);
			int x = (int)(p.x / L);
			//std::cout << p.x<<" "<<p.y<<" "<<L << " ";
			if (x < 0) { x = 0; }
			if (x > n - 1) { x = n - 1; }
			int y = (int)(p.y / L);
			if (y < 0) { y = 0; }
			if (y > n-1) { y = n-1; }
			Grids[x*n + y].input(p);
			//std::cout << x << " " << y << " ";
			//std::cout << x * n + y << std::endl;
		}

		inline void checkOrder() {
			//扫描一遍所有的格子，把有两个以上格子的点提出来
			for (int i = 0; i < Grids.size(); i++) {
				if (Grids[i].used < 2) {
					continue;
				}
				std::vector<std::vector<Point2>> p=Grids[i].points;
				int num_p = static_cast<int>(p.size());
				//std::cout << i<<" "<<num_p << std::endl;
				for (int j = 0; j < num_p; j++) {
					for (int k = j + 1; k < num_p; k++) {
						//std::cout << j << " " << k << std::endl;
						std::vector<Point2> points1 = p[j];
						std::vector<Point2> points2 = p[k];
						if (points1.size() > points2.size()) {
							std::swap(points1, points2);
						}
						int site1 = points1.front().father;
						int site2 = points2.front().father;

						int num = checkSitePair(site1, site2);
						if (num == -1) {
							SitePairPoints spp(site1, site2);
							for (int l = 0; l < points1.size(); l++) {
								spp.points1.push_back(points1[l]);
							}
							for (int l = 0; l < points2.size(); l++) {
								spp.points2.push_back(points2[l]);
							}
							sitePairsPoints.push_back(spp);
						}
						else {
							for (int l = 0; l < points1.size(); l++) {
								sitePairsPoints[num].points1.push_back(points1[l]);
							}
							for (int l = 0; l < points2.size(); l++) {
								sitePairsPoints[num].points2.push_back(points2[l]);
							}
						}
					}
				}
			}
			
			int sitePairsPointsSize = static_cast<int>(sitePairsPoints.size());
			for (int i = 0; i < sitePairsPointsSize; i++) {
				Entry* entry1 = &entries_[sitePairsPoints[i].site1];
				Entry* entry2 = &entries_[sitePairsPoints[i].site2];
				Site* site1 = entry1->site;
				Site* site2 = entry2->site;
				std::vector<Point2>* points1 = &sitePairsPoints[i].points1;
				std::vector<Point2>* points2 = &sitePairsPoints[i].points2; 
				
				typename Candidate::Vector candidates1(points1->size());
				int size = static_cast<int>(points1->size());
				for (int l = 0; l < size; l++) {
					Point& point = (*points1)[l];
					candidates1[l] = Candidate(&point, energy(point, site1), energy(point, site2));
				}
				std::make_heap(candidates1.begin(), candidates1.end());

				double minEnergy = -(candidates1.front().energySelf - candidates1.front().energyOther);

				typename Candidate::Vector candidates2(points2->size());
				size = static_cast<int>(points2->size());
				for (int l = 0; l < size; l++) {
					Point& point = (*points2)[l];
					double eSelf = energy(point, site2);
					double eOther = energy(point, site1);
					if (eSelf - eOther > minEnergy) {
						candidates2[l] = Candidate(&point, eSelf, eOther);
					}
				}
				std::make_heap(candidates2.begin(), candidates2.end());

				int maxSwaps = static_cast<int>(std::min(candidates1.size(), candidates2.size()));
				int swaps = 0;
				double deltaEnergy = 0;
				for (swaps = 0; swaps < maxSwaps; swaps++) {
					Candidate& candidate1 = candidates1.front();
					Candidate& candidate2 = candidates2.front();
					if (candidate1.energySelf - candidate1.energyOther + candidate2.energySelf - candidate2.energyOther <= 0) {
						break;
					}
					deltaEnergy += candidate1.energySelf - candidate1.energyOther + candidate2.energySelf - candidate2.energyOther;
					std::pop_heap(candidates1.begin(), candidates1.end() - swaps);
					std::pop_heap(candidates2.begin(), candidates2.end() - swaps);
				}
				SitePair sp(sitePairsPoints[i].site1, sitePairsPoints[i].site2, deltaEnergy);
				sitePairs.push(sp);
			}
		}

		inline int checkSitePair(const int site1, const int site2) {
			for (int i = 0; i < sitePairsPoints.size(); i++) {
				SitePairPoints spp = sitePairsPoints[i];
				if (spp.site1 == site1 && spp.site2 == site2) {
					return i;
				}
			}
			return -1;
		}

		inline bool myOptimize() {
			std::vector<bool> stability(entries_.size(), true);
			std::priority_queue<SitePair> sps;
			while (!sitePairs.empty()) {
				SitePair sp = sitePairs.top();
				sitePairs.pop();
				sps.push(sp);
				Entry* entry1 = &entries_[sp.site1];
				Entry* entry2 = &entries_[sp.site2];

				if (entry1->stable && entry2->stable ||
					metric_.distance(entry1->bounding.center, entry2->bounding.center) > entry1->bounding.radius + entry2->bounding.radius) {
					continue;
				}

				if (entry1->points.size() > entry2->points.size()) {
					std::swap(entry1, entry2);
				}

				Site* site1 = entry1->site;
				Site* site2 = entry2->site;
				//std::cout << site1->location.x << " " << site2->location.x << std::endl;
				std::vector<Point>* points1 = &entry1->points;
				std::vector<Point>* points2 = &entry2->points;

				double maxSquaredRadius = std::max(entry1->bounding.squaredRadius, entry2->bounding.squaredRadius);
				//printf("%d %d %f\n", i, j, maxSquaredRadius);

				typename Candidate::Vector candidates1(points1->size());
				int size = static_cast<int>(points1->size());
				int count = 0;
				for (int k = 0; k < size; ++k) {
					Point& point = (*points1)[k];
					if (metric_.distance_square(point, entry2->bounding.center) <= maxSquaredRadius) {
						candidates1[count++] = Candidate(&point, energy(point, site1), energy(point, site2));
					}
				}
				if (count == 0) {
					continue;
				}
				candidates1.resize(count);
				std::make_heap(candidates1.begin(), candidates1.end());
				//std::cout << "OK1"<<std::endl;

				double minEnergy = -(candidates1.front().energySelf - candidates1.front().energyOther);
				//std::cout << minEnergy << std::endl;
				typename Candidate::Vector candidates2(points2->size());
				size = static_cast<int>(points2->size());
				count = 0;
				for (int k = 0; k < size; ++k) {
					Point& point = (*points2)[k];
					if (metric_.distance_square(point, entry1->bounding.center) <= maxSquaredRadius) {
						double eSelf = energy(point, site2);
						double eOther = energy(point, site1);
						if (eSelf - eOther > minEnergy) {
							candidates2[count++] = Candidate(&point, eSelf, eOther);
						}
					}
				}
				if (count == 0) {
					continue;
				}
				candidates2.resize(count);
				std::make_heap(candidates2.begin(), candidates2.end());
				//std::cout << "OK2" << std::endl;

				int maxSwaps = static_cast<int>(std::min(candidates1.size(), candidates2.size()));
				int swaps;
				for (swaps = 0; swaps < maxSwaps; ++swaps) {
					Candidate& candidate1 = candidates1.front();
					Candidate& candidate2 = candidates2.front();
					if (candidate1.energySelf - candidate1.energyOther + candidate2.energySelf - candidate2.energyOther <= 0) {
						break;
					}
					std::swap(*candidate1.point, *candidate2.point);
					entry1->energy += candidate2.energyOther - candidate1.energySelf;
					entry2->energy += candidate1.energyOther - candidate2.energySelf;
					std::pop_heap(candidates1.begin(), candidates1.end() - swaps);
					std::pop_heap(candidates2.begin(), candidates2.end() - swaps);
				}

				if (swaps > 0) {
					num_swap += swaps;
					stability[sp.site1] = false;
					stability[sp.site2] = false;
					entry1->update(metric_);
					entry2->update(metric_);
				}
			}
			bool stable = true;
			for (int i = 0; i < entries_.size(); ++i) {
				entries_[i].stable = stability[i];
				stable &= stability[i];
			}
			std::swap(sps, sitePairs);
			return stable;
		}

	private:
	public:

		struct Bounding
		{

			Bounding()
				: radius(0), squaredRadius(0) {
			}

			Bounding(const Point& center, const double radius)
				: center(center), radius(radius), squaredRadius(radius * radius) {
			}

			void update(const Point& site, const typename std::vector<Point>& points, const Metric& metric) {
				center = site;
				squaredRadius = 0;
				int pointsSize = static_cast<int>(points.size());
				for (int i = 0; i < pointsSize; ++i) {
					squaredRadius = std::max(squaredRadius, metric.distance_square(center, points[i]));
					//printf("(%d,%d) (%d,%d) %d\n", center.x,center.y,points[i].x,points[i].y,metric.distance_square(center, points[i]));
				}
				radius = sqrt(squaredRadius);
			}

			Point   center;
			double  radius;
			double  squaredRadius;

		};

		struct Entry
		{

			typedef std::map<int, Entry*>	MapPtr;
			typedef	std::vector<Entry>    Vector;

			Entry()
				: site(NULL), stable(false) {
			}

			Entry(Site *const site)
				: bounding(site->location, 0), site(site), stable(false) {
			}

			void update(const Metric& metric) {
				stable = false;
				bounding.update(site->location, points, metric);
			}

			Bounding            bounding;
			std::vector<Point>  points;
			Site*               site;
			bool		        stable;
			double              energy;

		};

		struct Candidate
		{

			typedef std::vector<Candidate> Vector;

			Candidate()
				: point(NULL), energySelf(0), energyOther(0) {
			}

			Candidate(Point *const point, const double energySelf, const double energyOther)
				: point(point), energySelf(energySelf), energyOther(energyOther) {
			}

			inline bool operator<(const Candidate& candidate) const {
				return energySelf - energyOther < candidate.energySelf - candidate.energyOther;
			}

			Point*  point;
			double  energySelf;
			double  energyOther;

		};

		//The grid stores the points allocated to this grid 
		struct Grid {
			int num_x, num_y;
			int used;
			double L;
			std::vector<std::vector<Point2>> points;

			Grid(const int x, const int y, const double L) :num_x(x), num_y(y), L(L), used(0) {

			}

			void input(const Point2 p) {
				if (used != 0) {
					for (int i = 0; i < points.size(); i++) {
						if (points[i][0].father == p.father) {
							points[i].push_back(p);
							return;
						}
					}
					used++;
					std::vector<Point2> vp;
					vp.push_back(p);
					points.push_back(vp);
				}
				else {
					used++;
					std::vector<Point2> vp;
					vp.push_back(p);
					points.push_back(vp);
				}
			}
		};

		struct SitePair {
			int site1, site2;
			double deltaE;

			SitePair(const int site1, const int site2,const double deltaE) :site1(site1), site2(site2), deltaE(deltaE) {

			}

			inline bool operator<(const SitePair& sitePair) const {
				return deltaE < sitePair.deltaE;
			}
		};

		struct SitePairPoints {
			int site1, site2;
			std::vector<Point2> points1;
			std::vector<Point2> points2;

			SitePairPoints(const int site1, const int site2) :site1(site1), site2(site2) {

			}
		};

		inline double energy(const Point& point, const Site *const site) const {
			return metric_.distance_square(point, site->location);
		}

		int num_swap = 0;

		Metric                      metric_;
		typename Entry::Vector	    entries_;
		typename Entry::MapPtr      mapping_;
		typename std::vector<Site>  sites_;

		std::vector<Grid>				Grids;
		std::priority_queue<SitePair>	sitePairs;//site
		std::vector<SitePairPoints>		sitePairsPoints;
	};

}

#endif	// CCVT_OPTIMIZER_H