#ifndef CCVT_SITE_H
#define CCVT_SITE_H

#include <list>
#include <vector>

namespace ccvt {

	template<class Point>
	struct Site
	{
		typedef std::list<Site> List;
		typedef std::vector<Site> Vector;
		typedef std::vector<Site*> VectorPtr;

		Site() :id(-1), capacity(0) {
		}

		Site(const int id, const int capacity, const Point& location) :id(id), capacity(capacity), location(location) {
		}

		int   id;//վ����
		int	  capacity;//վ������
		Point location;//վ��λ��
	};
}

#endif // !1
