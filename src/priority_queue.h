#pragma once

#include <map>
#include <queue>

template <typename T, typename Compare = std::less<T>>
class PriorityQueue {

    struct I {
        size_t index;
        size_t date;
    };

public:
    PriorityQueue(const Compare& cmp) {
        auto pair_cmp = [cmp](std::pair<T, I> a, std::pair<T, I> b) {
            return cmp(a.first, b.first);
        };
        q_ = std::priority_queue<std::pair<T, I>, std::vector<std::pair<T, I>>, std::function<bool(std::pair<T, I>, std::pair<T, I>)>>(pair_cmp);
    }

    bool empty() {
        while (!q_.empty()) {
            std::pair<T, I> x = q_.top();
            if (x.second.date == latest_dates_[x.second.index]) {
                return false;
            }
            q_.pop();
        }
        return true;
    }

    void push_or_update(T value, size_t index) {
        q_.push({value, {index, ++latest_dates_[index]}});
    }

    T pop() {
        while (!q_.empty()) {
            std::pair<T, I> x = q_.top();
            q_.pop();
            if (x.second.date == latest_dates_[x.second.index]) {
                return x.first;
            }
        }
        throw std::logic_error("Pop on empty priority queue");
    }

private:
    std::priority_queue<std::pair<T, I>, std::vector<std::pair<T, I>>, std::function<bool(std::pair<T, I>, std::pair<T, I>)>> q_;

    std::map<size_t, size_t> latest_dates_;
};
