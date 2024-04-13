#pragma once
#include <vector>
#include <span>
#include <utility>

namespace dashing2 {

// Source: https://gist.github.com/Tmplt/60ec8e3923bf6c09ffff655fbc3c3485
// Modified to use signed integers.

template<typename V>
static inline size_t edit_distance(const std::span<V> s1, const std::span<V> s2)
{
    int64_t len1 = s1.size();
    int64_t len2 = s2.size();

    int64_t i, j;

    /* Base case: a string is empty. */
    if (len1 == 0) return len2;
    if (len2 == 0) return len1;

    /*
     * Initialize rows[i].first (the previous row of distances).
     * This first row represents the edit distance againts an
     * empty string, so the distance is just the numbers of
     * characters to delete from the non-empty string.
     */
    std::vector<std::pair<int64_t, int64_t>> rows(len2 + 1);
    for (i = 0; i < static_cast<int64_t>(rows.size()); i++)
        rows[i].first = i;

    for (i = 0; i < len1; i++) {
        /* Edit distance is delete (i+1) chars from non-empty to empty. */
        rows[0].second = i + 1;

        /* Fill in the rest of the row. */
        for (j = 0; j < len2; j++) {
            auto cost = (s1[i] == s2[j]) ? 0 : 1;

            rows[j + 1].second = std::min({rows[j].second + 1,
                                           rows[j + 1].first + 1,
                                           rows[j].first + cost});
        }

        /*
         * Copy the current row to the previous
         * one in preparation for next iteration.
         */
        for (auto &row : rows)
            row.first = row.second;
    }

    return rows[len2].second;
}

}
