/* Copyright (c) 1997-2018
   Ewgenij Gawrilow, Michael Joswig (Technische Universitaet Berlin, Germany)
   http://www.polymake.org

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version: http://www.gnu.org/licenses/gpl.txt.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
--------------------------------------------------------------------------------
*/

#ifndef POLYMAKE_RATIONAL_FUNCTION_H
#define POLYMAKE_RATIONAL_FUNCTION_H

#include "polymake/Polynomial.h"

namespace pm {

template <typename T, typename Coefficient, typename Exponent>
struct is_unipolynomial_type {
   static const bool value=is_derived_from<T, UniPolynomial<Coefficient, Exponent> >::value;
};

template <typename Coefficient=Rational, typename Exponent=int>
class RationalFunction {
public:
   typedef UniPolynomial<Coefficient, Exponent> polynomial_type;
   typedef typename polynomial_type::term_hash term_hash;
   typedef Coefficient coefficient_type;
   typedef Exponent exponent_type;

   template <typename T>
   using fits_as_coefficient = typename polynomial_type::template fits_as_coefficient<T>;

   template <typename T>
   struct fits_as_particle
      : bool_constant<fits_as_coefficient<T>::value || is_unipolynomial_type<T, Coefficient, Exponent>::value> {};

   /// Construct a zero value.
   RationalFunction()
      : num()
      , den(one_value<coefficient_type>()) {}

   /// Construct a value with denominator equal to 1.
   template <typename T, typename enabled=typename std::enable_if<fits_as_particle<T>::value>::type>
   explicit RationalFunction(const T& c)
      : num(c)
      , den(one_value<coefficient_type>()) {}

   /// Construct a value with denominator equal to a constant scalar.
   template <typename T1, typename T2>
   RationalFunction(const T1& c1, const T2& c2,
                    typename std::enable_if<fits_as_particle<T1>::value && fits_as_coefficient<T2>::value, void**>::type=nullptr)
      : num(c1)
      , den(one_value<coefficient_type>())
   {
      num /= c2;
   }

   /// Construct a value with numerator equal to a constant scalar.
   template <typename T1, typename T2>
   RationalFunction(const T1& c1, const T2& p2,
                    typename std::enable_if<fits_as_coefficient<T1>::value &&
                                            is_unipolynomial_type<T2, Coefficient, Exponent>::value, void**>::type=nullptr)
      : num(c1)
      , den(p2)
   {
      if (is_zero(p2)) throw GMP::ZeroDivide();
      normalize_lc();
   }

   /// Construct a value with given numerator and denominator.
   template <typename T1, typename T2>
   RationalFunction(const T1& p1, const T2& p2,
                    typename std::enable_if<is_unipolynomial_type<T1, Coefficient, Exponent>::value &&
                                            is_unipolynomial_type<T2, Coefficient, Exponent>::value, void**>::type=nullptr)
   {
      if (is_zero(p2)) throw GMP::ZeroDivide();
      simplify(p1, p2);
      normalize_lc();
   }

   typedef polynomial_type numerator_type;
   typedef polynomial_type denominator_type;

   friend
   const polynomial_type& numerator(const RationalFunction& me) { return me.num; }

   friend
   const polynomial_type& denominator(const RationalFunction& me) { return me.den; }

   void swap(RationalFunction& other)
   {
      num.swap(other.num);
      den.swap(other.den);
   }

   RationalFunction& negate()
   {
      num.negate();
      return *this;
   }

   friend
   RationalFunction operator- (const RationalFunction& me)
   {
      return RationalFunction(-me.num, me.den, std::true_type());
   }

   template <typename T>
   typename std::enable_if<fits_as_particle<T>::value, RationalFunction&>::type
   operator+= (const T& c)
   {
      if (__builtin_expect(!is_zero(c), 1))
         num += den * c;
      return *this;
   }

   RationalFunction& operator+= (const RationalFunction& r)
   {
      if (__builtin_expect(!r.num.trivial(), 1)) {
         ExtGCD<polynomial_type> x = ext_gcd(den, r.den, false);
         x.p = x.k1 * x.k2; // x.p is used as dummy variable
         den.swap(x.p);
         x.k1 *= r.num;  x.k1 += num * x.k2;
         if( !is_one(x.g) ){
            x = ext_gcd(x.k1, x.g);
            x.k2 *= den;
            den.swap(x.k2);
         }
         num.swap(x.k1);
         normalize_lc();
      }
      return *this;
   }

   template <typename T>
   typename std::enable_if<fits_as_particle<T>::value, RationalFunction&>::type
   operator-= (const T& r)
   {
      if (__builtin_expect(!is_zero(r), 1))
         num -= den * r;
      return *this;
   }

   RationalFunction& operator-= (const RationalFunction& r)
   {
      if (__builtin_expect(!r.num.trivial(), 1)) {
         ExtGCD<polynomial_type> x = ext_gcd(den, r.den, false);
         x.p = x.k1 * x.k2; // x.p is used as dummy variable
         den.swap(x.p);
         x.k1 *= r.num;  x.k1.negate();  x.k1 += num * x.k2;
         if( !is_one(x.g) ){
            x = ext_gcd(x.k1, x.g);
            x.k2 *= den;
            den.swap(x.k2);
         }
         num.swap(x.k1);
         normalize_lc();
      }
      return *this;
   }

   template <typename T> friend
   typename std::enable_if<fits_as_particle<T>::value, RationalFunction>::type
   operator+ (const RationalFunction& l, const T& r)
   {
      if (__builtin_expect(!is_zero(r), 1))
         return RationalFunction(l.num + l.den * r, l.den, std::true_type());
      else
         return l;
   }

   template <typename T> friend
   typename std::enable_if<fits_as_particle<T>::value, RationalFunction>::type
   operator+ (const T& l, const RationalFunction& r)
   {
      return r+l;
   }

   template <typename T> friend
   typename std::enable_if<fits_as_particle<T>::value, RationalFunction>::type
   operator- (const RationalFunction& l, const T& r)
   {
      if (__builtin_expect(!is_zero(r), 1))
         return RationalFunction(l.num - l.den * r, l.den, std::true_type());
      else
         return l;
   }

   template <typename T> friend
   typename std::enable_if<fits_as_particle<T>::value, RationalFunction>::type
   operator- (const T& l, const RationalFunction& r)
   {
      if (__builtin_expect(!is_zero(l), 1))
         return RationalFunction(r.den * l - r.num, r.den, std::true_type());
      else
         return -r;
   }

private:
   RationalFunction& normalize_after_addition(ExtGCD<polynomial_type>& x)
   {
      if (!is_one(x.g)) {
         x = ext_gcd(num, x.g);
         x.k2 *= den;
         den.swap(x.k2);
         num.swap(x.k1);
      }
      normalize_lc();
      return *this;
   }

public:
   friend
   RationalFunction operator+ (const RationalFunction& l, const RationalFunction& r)
   {
      if (__builtin_expect(l.num.trivial(), 0)) {
         return r;
      } else if (__builtin_expect(r.num.trivial(), 0)) {
         return l;
      } else {
         ExtGCD<polynomial_type> x = ext_gcd(l.den, r.den, false);
         return RationalFunction(l.num * x.k2 + r.num * x.k1, x.k1 * x.k2, std::true_type()).normalize_after_addition(x);
      }
   }

   friend
   RationalFunction operator- (const RationalFunction& l, const RationalFunction& r)
   {
      if (__builtin_expect(l.num.trivial(), 0)) {
         return -r;
      } else if (__builtin_expect(r.num.trivial(), 0)) {
         return l;
      } else {
         ExtGCD<polynomial_type> x = ext_gcd(l.den, r.den, false);
         return RationalFunction(l.num * x.k2 - r.num * x.k1, x.k1 * x.k2, std::true_type()).normalize_after_addition(x);
      }
   }

   template <typename T>
   typename std::enable_if<fits_as_coefficient<T>::value, RationalFunction>::type&
   operator*= (const T& r)
   {
      num *= r;
      return *this;
   }

   template <typename T>
   typename std::enable_if<fits_as_coefficient<T>::value, RationalFunction>::type&
   operator/= (const T& r)
   {
      num /= r;
      return *this;
   }

   template <typename T> friend
   typename std::enable_if<fits_as_coefficient<T>::value, RationalFunction>::type
   operator* (const RationalFunction& l, const T& r)
   {
      if (__builtin_expect(!is_zero(r), 1))
         return RationalFunction(l.num * r, l.den, std::true_type());
      else
         return RationalFunction();
   }

   template <typename T> friend
   typename std::enable_if<fits_as_coefficient<T>::value, RationalFunction>::type
   operator* (const T& l, const RationalFunction& r)
   {
      if (__builtin_expect(!is_zero(l), 1))
         return RationalFunction(l * r.num, r.den, std::true_type());
      else
         return RationalFunction();
   }

   template <typename T> friend
   typename std::enable_if<fits_as_coefficient<T>::value, RationalFunction>::type
   operator/ (const RationalFunction& l, const T& r)
   {
      return RationalFunction(l.num / r, l.den, std::true_type());
   }

   template <typename T> friend
   typename std::enable_if<fits_as_coefficient<T>::value, RationalFunction>::type
   operator/ (const T& l, const RationalFunction& r)
   {
      if (__builtin_expect(r.num.trivial(), 0)) {
         throw GMP::ZeroDivide();
      } else if (__builtin_expect(is_zero(l), 0)) {
         return RationalFunction();
      } else {
         return RationalFunction(l * r.den, r.num, std::false_type());
      }
   }

   template <typename T> friend
   typename std::enable_if<is_unipolynomial_type<T, Coefficient, Exponent>::value, RationalFunction>::type
   operator* (const RationalFunction& l, const T& r)
   {
      if (__builtin_expect(is_zero(r), 1)) {
         return RationalFunction();
      } else if (__builtin_expect(l.num.trivial(), 0)) {
         return l;
      } else {
         const ExtGCD<polynomial_type> x = ext_gcd(l.den, r, false);
         return RationalFunction(l.num * x.k2, x.k1, std::true_type());
      }
   }

   template <typename T> friend
   typename std::enable_if<is_unipolynomial_type<T, Coefficient, Exponent>::value, RationalFunction>::type
   operator* (const T& l, const RationalFunction& r)
   {
      return r * l;
   }

   friend
   RationalFunction operator* (const RationalFunction& l, const RationalFunction& r)
   {
      if (__builtin_expect(l.num.trivial(), 0)) {
         return l;
      } else if (__builtin_expect(r.num.trivial(), 0)) {
         return r;
      } else if (l.den==r.den || l.num==r.num) {
         return RationalFunction(l.num * r.num, l.den * r.den, std::true_type());
      } else {
         const ExtGCD<polynomial_type> x = ext_gcd(l.num, r.den, false),
                                       y = ext_gcd(l.den, r.num, false);
         return RationalFunction(x.k1 * y.k2, y.k1 * x.k2, std::false_type());
      }
   }

   template <typename T> friend
   typename std::enable_if<is_unipolynomial_type<T, Coefficient, Exponent>::value, RationalFunction>::type
   operator/ (const RationalFunction& l, const T& r)
   {
      if (__builtin_expect(is_zero(r), 1)) {
         throw GMP::ZeroDivide();
      } else if (__builtin_expect(l.num.trivial(), 0)) {
         return l;
      } else {
         const ExtGCD<polynomial_type> x = ext_gcd(l.num, r, false);
         return RationalFunction(x.k1, l.den * x.k2, std::false_type());
      }
   }

   template <typename T> friend
   typename std::enable_if<is_unipolynomial_type<T, Coefficient, Exponent>::value, RationalFunction>::type
   operator/ (const T& l, const RationalFunction& r)
   {
      if (__builtin_expect(r.num.trivial(), 1)) {
         throw GMP::ZeroDivide();
      } else if (__builtin_expect(is_zero(l), 0)) {
         return RationalFunction();
      } else {
         const ExtGCD<polynomial_type> x = ext_gcd(r.num, l, false);
         return RationalFunction(r.den * x.k2, x.k1, std::false_type());
      }
   }

   friend
   RationalFunction operator/ (const RationalFunction& l, const RationalFunction& r)
   {
      if (__builtin_expect(r.num.trivial(), 0)) {
         throw GMP::ZeroDivide();
      } else if (__builtin_expect(l.num.trivial(), 0)) {
         return l;
      } else if (l.den==r.num || l.num==r.den) {
         return RationalFunction(l.num * r.den, l.den * r.num, std::true_type()); 
      } else {
         const ExtGCD<polynomial_type> x = ext_gcd(l.num, r.num, false),
                                       y = ext_gcd(l.den, r.den, false);
         return RationalFunction(x.k1 * y.k2, y.k1 * x.k2, std::false_type());
      }
   }

   template <typename T>
   typename std::enable_if<is_unipolynomial_type<T, Coefficient, Exponent>::value, RationalFunction>::type&
   operator*= (const T& r)
   {
      *this = (*this) * r;
      return *this;
   }

   template <typename T>
   typename std::enable_if<is_unipolynomial_type<T, Coefficient, Exponent>::value, RationalFunction>::type&
   operator/= (const T& r)
   {
      *this = (*this) / r;
      return *this;
   }

   RationalFunction& operator*= (const RationalFunction& r)
   {
      *this = (*this) * r;
      return *this;
   }

   RationalFunction& operator/= (const RationalFunction& r)
   {
      *this = (*this) / r;
      return *this;
   }

   friend
   bool operator== (const RationalFunction& l, const RationalFunction& r)
   {
      return l.num == r.num && l.den == r.den;
   }

   template <typename T> friend
   typename std::enable_if<fits_as_particle<T>::value, bool>::type
   operator== (const RationalFunction& l, const T& r)
   {
      return l.den.is_one() && l.num == r;
   }

   template <typename T> friend
   typename std::enable_if<fits_as_particle<T>::value, bool>::type
   operator== (const T& l, const RationalFunction& r)
   {
      return r == l;
   }

   friend
   bool operator!= (const RationalFunction& l, const RationalFunction& r)
   {
      return !(l == r);
   }

   template <typename T> friend
   typename std::enable_if<fits_as_particle<T>::type, bool>::type
   operator!= (const RationalFunction& l, const T& r)
   {
      return !(l == r);
   }

   template <typename T> friend
   typename std::enable_if<fits_as_particle<T>::value, bool>::type
   operator!= (const T& l, const RationalFunction& r)
   {
      return !(r == l);
   }

   static
   const Array<std::string>& get_var_names()
   {
      return polynomial_type::get_var_names();
   }

   static
   void set_var_names(const Array<std::string>& names)
   {
      polynomial_type::set_var_names(names);
   }

   static
   void reset_var_names()
   {
      polynomial_type::reset_var_names();
   }

   static
   void swap_var_names(PolynomialVarNames& other_names)
   {
      polynomial_type::swap_var_names(other_names);
   }

   template <typename Output> friend
   Output& operator<< (GenericOutput<Output>& out, const RationalFunction& rf)
   {
      out.top() << '(' << rf.num << ")/(" << rf.den << ')';
      return out.top();
   }

   explicit operator const polynomial_type& () const
   {
      if (!den.is_one())
         throw std::runtime_error("Denominator is not one; cannot convert to a polynomial");
      return num;
   }

   size_t get_hash() const noexcept
   {
      size_t h=num.get_hash();
      hash_combine(h, den.get_hash());
      return h;
   }

protected:
   polynomial_type num, den;

   /// internally used constructor, when the operands are known to be mutually prime and normalized
   template <bool lc_normalized>
   RationalFunction(const polynomial_type& num_arg, const polynomial_type& den_arg, bool_constant<lc_normalized>) :
      num(num_arg),
      den(den_arg)
   {
      if (!lc_normalized) normalize_lc();
   }

   void normalize_lc()
   {
      if (num.trivial()) {
         den=polynomial_type(one_value<coefficient_type>());
      } else {
         const Coefficient den_lc=den.lc();
         if (!is_one(den_lc)) {
            num /= den_lc;
            den /= den_lc;
         }
      }
   }

   // constructor helpers

   void simplify(const Coefficient& c1, const Exponent& e1,
                 const Coefficient& c2, const Exponent& e2)
   {
      if (e1 < e2)
      {
         // x^e1 / x^e2  ==  1 / x^(e2-e1)
         num=polynomial_type(c1);
         den=polynomial_type(e2-e1, c2);
      }
      else
      {
         // x^e1 / x^e2  ==  x^(e1-e2)
         num=polynomial_type(e1-e2, c1);
         den=polynomial_type(c2);
      }
   }

   void simplify(const polynomial_type& p1,
                 const Coefficient& c2, const Exponent& e2)
   {
      const Exponent e1=p1.lower_deg();
      if (e1 < e2)
      {
         // r(x)*x^e1 / x^e2  ==  r(x) / x^(e2-e1)
         if ( !is_zero(e1) ) {
            div_exact(p1, polynomial_type(e1, one_value<coefficient_type>())).swap(num);
         } else {
            num=p1;
         }
         den=polynomial_type(polynomial_type(e2-e1, one_value<coefficient_type>()), c2);
      }
      else
      {
         // r(x)*x^e1 / x^e2  ==  r(x)*x^(e1-e2)
         div_exact(p1, polynomial_type(e2, one_value<coefficient_type>())).swap(num);
         den=polynomial_type(c2);
      }
   }

   void simplify(const polynomial_type& p1, const polynomial_type& p2)
   {
      ExtGCD<polynomial_type> x=ext_gcd(p1, p2, false);
      num.swap(x.k1);
      den.swap(x.k2);
   }

   template <typename> friend struct spec_object_traits;

#if POLYMAKE_DEBUG
public:
   void dump() const __attribute__((used)) { cerr << *this << std::flush; }
#endif
};

template <typename Coefficient, typename Exponent>
struct spec_object_traits< Serialized< RationalFunction<Coefficient, Exponent> > > :
   spec_object_traits<is_composite> {

   typedef RationalFunction<Coefficient, Exponent> masquerade_for;

   typedef cons<typename RationalFunction<Coefficient, Exponent>::term_hash,
                typename RationalFunction<Coefficient, Exponent>::term_hash> elements;

   template <typename Me, typename Visitor>
   static void visit_elements(Me& me, Visitor& v)
   {
      v << me.num.get_mutable_terms() << me.den.impl_ptr->get_mutable_terms();
   }

};

template <typename Coefficient, typename Exponent>
struct choose_generic_object_traits< RationalFunction<Coefficient, Exponent>, false, false > :
      spec_object_traits< RationalFunction<Coefficient, Exponent> > {
   typedef void generic_type;
   typedef is_polynomial generic_tag;
   typedef RationalFunction<Coefficient, Exponent> persistent_type;

   static
   bool is_zero(const persistent_type& p)
   {
      return numerator(p).trivial();
   }

   static
   bool is_one(const persistent_type& p)
   {
      return numerator(p).is_one() && denominator(p).is_one();
   }

   static
   const persistent_type& zero()
   {
      static const persistent_type x=persistent_type();
      return x;
   }

   static
   const persistent_type& one()
   {
      static const persistent_type x(1);
      return x;
   }
};

template <typename Coefficient, typename Exponent>
struct algebraic_traits< RationalFunction<Coefficient, Exponent> > {
   typedef RationalFunction<typename algebraic_traits<Coefficient>::field_type, Exponent> field_type;
};

namespace polynomial_impl {

template <typename Coefficient, typename Exponent>
struct nesting_level< RationalFunction<Coefficient, Exponent> >
  : int_constant<nesting_level<Coefficient>::value+1> {};

}

template <typename Coefficient, typename Exponent, typename T, typename TModel>
struct isomorphic_types_impl<RationalFunction<Coefficient, Exponent>, T,
                             typename std::enable_if<RationalFunction<Coefficient, Exponent>::template fits_as_particle<T>::value, is_polynomial>::type,
                             TModel>
   : std::false_type {
   typedef cons<is_polynomial, is_scalar> discriminant;
};

template <typename Coefficient, typename Exponent, typename T, typename TModel>
struct isomorphic_types_impl<T, RationalFunction<Coefficient, Exponent>, TModel,
                             typename std::enable_if<RationalFunction<Coefficient, Exponent>::template fits_as_particle<T>::value, is_polynomial>::type>
   : std::false_type {
   typedef cons<is_scalar, is_polynomial> discriminant;
};

template <typename Coefficient, typename Exponent>
struct isomorphic_types_impl<RationalFunction<Coefficient, Exponent>, RationalFunction<Coefficient, Exponent>, is_polynomial, is_polynomial>
   : std::true_type {
   typedef cons<is_polynomial, is_polynomial> discriminant;
};


template <typename Coefficient, typename Exponent, typename T> inline
typename std::enable_if<is_unipolynomial_type<T, Coefficient, Exponent>::value, RationalFunction<Coefficient, Exponent>>::type
operator/ (const UniPolynomial<Coefficient, Exponent>& p1, const T& p2)
{
   return RationalFunction<Coefficient, Exponent>(p1, p2);
}

template <typename Coefficient, typename Exponent, typename T> inline
typename std::enable_if<RationalFunction<Coefficient, Exponent>::template fits_as_coefficient<T>::value,
                        RationalFunction<Coefficient, Exponent>>::type
operator/ (const T& p1, const UniPolynomial<Coefficient, Exponent>& p2)
{
   return RationalFunction<Coefficient, Exponent>(p1, p2);
}

namespace operations {

// operations neg, add, sub, mul defined in Polynomial.h will instantiate correctly for RationalFunction too.

template <typename LeftRef, class RightRef>
struct div_impl<LeftRef, RightRef,
                typename std::enable_if<is_instance_of<typename deref<LeftRef>::type, RationalFunction>::value, cons<is_polynomial, is_polynomial>>::type> {
   typedef LeftRef first_argument_type;
   typedef RightRef second_argument_type;
   typedef typename deref<LeftRef>::type result_type;

   result_type operator() (typename function_argument<LeftRef>::const_type l,
                           typename function_argument<RightRef>::const_type r) const
   {
      return l/r;
   }

   void assign(typename lvalue_arg<LeftRef>::type l, typename function_argument<RightRef>::const_type r) const
   {
      l/=r;
   }
};

template <typename LeftRef, class RightRef>
struct div_impl<LeftRef, RightRef,
                typename std::enable_if<is_instance_of<typename deref<RightRef>::type, RationalFunction>::value, cons<is_scalar, is_polynomial>>::type> {
   typedef LeftRef first_argument_type;
   typedef RightRef second_argument_type;
   typedef typename deref<RightRef>::type result_type;

   result_type operator() (typename function_argument<LeftRef>::const_type l,
                           typename function_argument<RightRef>::const_type r) const
   {
      return l/r;
   }
};

template <typename LeftRef, class RightRef>
struct div_impl<LeftRef, RightRef,
                typename std::enable_if<is_instance_of<typename deref<LeftRef>::type, RationalFunction>::value, cons<is_polynomial, is_scalar>>::type> {
   typedef LeftRef first_argument_type;
   typedef RightRef second_argument_type;
   typedef typename deref<LeftRef>::type result_type;

   result_type operator() (typename function_argument<LeftRef>::const_type l,
                           typename function_argument<RightRef>::const_type r) const
   {
      return l/r;
   }

   void assign(typename lvalue_arg<LeftRef>::type l, typename function_argument<RightRef>::const_type r) const
   {
      l/=r;
   }
};

} // end namespace operations 

} // end namespace pm

namespace polymake {
   using pm::RationalFunction;
}

namespace std {
   template <typename Coefficient, typename Exponent>
   void swap(pm::RationalFunction<Coefficient,Exponent>& x1, pm::RationalFunction<Coefficient,Exponent>& x2) { x1.swap(x2); }
}

#endif // POLYMAKE_RATIONAL_FUNCTION_H

// Local Variables:
// mode:C++
// c-basic-offset:3
// indent-tabs-mode:nil
// End:
