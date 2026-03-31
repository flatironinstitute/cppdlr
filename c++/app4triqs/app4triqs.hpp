/**
 * @file
 * @brief Provides various functions and classes for the **app4triqs** application.
 */

#pragma once
#include <h5/h5.hpp>

namespace app4triqs {

  /**
   * @addtogroup mygroup
   * @{
   */

  /**
   * A very useful and important class
   *
   * @note A Useful note
   */
  class toto {

    int i = 0;

    public:
    toto() = default;

    /**
     * Construct from integer
     *
     * @param i_ a scalar \f$ G(\tau) \f$
     */
    explicit toto(int i_) : i(i_) {}

    ~toto() = default;

    // Copy/Move construction
    toto(toto const &) = default;
    toto(toto &&)      = default;

    /// Copy/Move assignment
    toto &operator=(toto const &) = default;
    toto &operator=(toto &&)      = default;

    /// Simple accessor
    [[nodiscard]] int get_i() const { return i; }

    /** 
     * A simple function with \f$ G(\tau) \f$
     *
     * @param u Nothing useful
     */
    int f(int u) { return u; }

    /// Arithmetic operations
    toto operator+(toto const &b) const;
    toto &operator+=(toto const &b);

    /// Comparison
    bool operator==(toto const &b) const;

    /// HDF5
    static std::string hdf5_format() { return "Toto"; }

    friend void h5_write(h5::group grp, std::string subgroup_name, toto const &m);
    friend void h5_read(h5::group grp, std::string subgroup_name, toto &m);

    /// Serialization
    void serialize(auto &ar) const { ar & i; }
    void deserialize(auto &ar) { ar & i; }
  };

  /**
   * @brief Chain digits of two integers
   *
   * @details A set of functions that implement chaining
   *
   * Do I really need to explain more ? 
   *
   * @param i The first integer
   * @param j The second integer
   * @return An integer containing the digits of both i and j
   */
  int chain(int i, int j);

  /** @} */

} // namespace app4triqs
