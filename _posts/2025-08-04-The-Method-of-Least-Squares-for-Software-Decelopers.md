---
title: The Method of Least Squares for Software Developers
layout: post
permalink: /the-method-of-least-squares-for-software-developers/
---

* Toc placeholder
{:toc}

## Introduction

I am a seasoned C++ developer and in many of my jobs I also interviewed job candidates. I normally interview for C++ but when the role involves graphics, I also interview for Linear Algebra. One of the Linear Algebra problems I give to candidates is the following one:

> You've got a plane defined by 3 arbitrary (not colinear) points on that plane. Let's call them $ \mathbf{P_1} $, $ \mathbf{P_2} $ and $ \mathbf{P_3} $. You've also got a 3rd point $ \mathbf{Q} $ that's generally not on the plane. Your goal is to find the point $ \mathbf{Q'} $ on that plane that's closest (in the Euclidean sense) to point $ \mathbf{Q} $. In other words, I am asking you to drop a perpendicular from $ \mathbf{Q} $ to the plane and find the intersection point of that perpendicular with the plane.

What I like about this problem is that there are many ways to solve it - I know at least 4. About every other candidate manages to solve it, yet over the years, no candidate has ever used my favourite approach. That's a shame, because that approach works for projecting to any-dimensional linear surface in any-dimensional space. So, you can use it to project not just to planes but also to lines in both 2D and 3D.

Given that the method in question seems to be little known among the software developers (yet very well known by non-software engineers), I though it would be a good idea to write a post about it.

## The Method of Least Squares

Or alternatively, you may know it as *Projecting Onto a Subspace*. If you haven't heard of either of those names, yet you know some basic Linear Algebra, this article is for you.

So, what does the method of least squares do? Consider an overdetermined linear system $ \mathbf{Ax} = \mathbf{b} $. By overdetermined, I mean the matrix $ \mathbf{A} $ is tall - it has more rows than columns. Generally, such linear systems don't have exact solutions. The best we can do is to find an approximate solution that minimizes the error. That error is called the residual, and that residual is equal to $ \mathbf{b} - \mathbf{Ax} $. The residual is a vector, so what do we mean by minimizing it? Well, we minimize its Euclidean norm:

$$
argmin_{\mathbf{x}} \| \mathbf{b} - \mathbf{Ax} \|
$$

Given that the Euclidean norm is non-negative, minimizing it is equivalent to minimizing its square:

$$
argmin_{\mathbf{x}} \| \mathbf{b} - \mathbf{Ax} \|^2
$$

Minimizing the squared norm is much easier from the Calculus point of view - the derivatives will be much nicer. So, we'll be minimizing the squared norm of the residual. Hence the *Method of Least Squares* name.

One other interpretation of this method is the following: it finds a point (techically a vector) in the column space of matrix $ \mathbf{A} $ that's closest in the Euclidean sense to point $ \mathbf{b} $. What's a column space of a matrix? It's a set of all possible linear combinations of columns of that matrix. In general, a set of all possible linear combinations of a given set of vectors is called a *vector space*.

How does a column space of a 3x2 (3 rows, 2 columns) matrix looks like? Well, it looks like a plane in 3D (provided the columns are linearly independent). Isn't that what we are after? We are trying to find a point on a plane that's closest in the Euclidean sense to a given point. Sounds like exactly what we want!

There is one tiny complication though: not every plane in 3D is a vector space, but only those that pass through the origin. Why is that? Well, a linear combination with all coefficients set to zero will produce a zero vector, so the zero vector (the origin) has to be a part of any vector space. This complication is easy to overcome though: we'll just shift our coordinate system to make one of the points on the plane (say $ \mathbf{P_1} $) to be our new origin. Then we solve the problem in that shifted coordinate system and then shift the result back into the original coordinate system.

Applying the method of least squares is easy: you just take the original linear system $ \mathbf{Ax} = \mathbf{b} $ and left-multiply both sides by $ \mathbf{A}^\top $. That gives you the so called *Normal Equation*:

$$
\mathbf{A}^\top \mathbf{Ax} = \mathbf{A}^\top \mathbf{b}
$$

When applied to our problem, we have:

$$
\begin{align}
\mathbf{A} &= \left[
\begin{array}{c|c}
\mathbf{P_2} - \mathbf{P_1} & \mathbf{P_3} - \mathbf{P_1}
\end{array} \notag \\
\right] \\
\mathbf{b} &= \mathbf{Q} - \mathbf{P_1} \notag
\end{align}
$$

Subtracting $ \mathbf{P_1} $ just shifts the coordinate system so that point $ \mathbf{P_1} $ is the new origin and that origin is on our plane.

Then we solve the normal equation for $ \mathbf{x} $. Once we've got $ \mathbf{x} $, we can get the solution in the original coordinate space by evaluating $ \mathbf{Ax} + \mathbf{P_1} $. Adding $ \mathbf{P_1} $ shifts us back into the original coordinate system.

Putting everything together, the solution to our original problem is going to be:

$$
\begin{align}
\mathbf{A} &= \left[
\begin{array}{c|c}
\mathbf{P_2} - \mathbf{P_1} & \mathbf{P_3} - \mathbf{P_1}
\end{array}
\right] \notag \\
\mathbf{b} &= \mathbf{Q} - \mathbf{P_1} \notag \\
\mathbf{x} &= (\mathbf{A}^\top \mathbf{A})^{-1}\mathbf{A}^\top\mathbf{b} \notag \\
\mathbf{Q'} &= \mathbf{Ax} + \mathbf{P_1} \notag
\end{align}
$$

When projecting onto a plane, $ \mathbf{A}^\top\mathbf{A} $ is going to be a 2x2 matrix. When projecting onto a line (whether in 3D or in 2D), it would be a 1x1 matrix, which makes inverting it really easy - just take the reciprocal of its only element.

## Show Me The Code!

Here is a function to project a point onto a plane in C++ using the [glm](https://github.com/g-truc/glm) library:

```c++
glm::vec3 projectToPlane(glm::vec3 p1, glm::vec3 p2, glm::vec3 p3, glm::vec3 q)
{
    // Note that in glm, a 2x3 matrix means a 2-column, 3-row, matrix.
    // In Linear Algebra, we would have called such a matrix a 3x2 one.
    const glm::mat2x3 A(p2 - p1, p3 - p1);

    const auto b = q - p1;
    const auto At(glm::transpose(A));
    const auto invAtA(glm::inverse(At * A));
    const auto x = invAtA * At * b;

    return p1 + A * x;
}
```

You can now project points onto a plane or onto a line (with minimal modifications to the code above). However, understanding how the normal equation is derived will allow you to solve a wider class of problems, one of which is given towards the end of this article. So, let's do the derivation!

## Deriving the Normal Equation

So, we want to project something onto a vector space. Let's recall that a vector space is a set of all possible linear combinations of a given set of vectors. Those vectors are called the basis vectors and in our case they are the columns of matrix $ \mathbf{A} $. Now suppose we also have matrix $ \mathbf{B} $, whose column space is an *orthogonal complimentary* one with respect to the column space of $ \mathbf{A} $. What does that mean exactly? The *orthogonal* part means that every column of $ \mathbf{A} $ is orthogonal to every column of $ \mathbf{B} $. The *complimentary* part means that when columns of $ \mathbf{A} $ and $ \mathbf{B} $ are put together, the vector space they produce covers the whole space (in this case, the whole 3D space). That is, every point in 3D shall be representable as a linear combination of columns of $ \mathbf{A} $ and $ \mathbf{B} $.

When projecting to a plane in 3D, we can actually build that matrix $ \mathbf{B} $ pretty easily. We just put the cross product of $ (\mathbf{P_2} - \mathbf{P_1}) $ and $ (\mathbf{P_3} - \mathbf{P_1}) $ as its only column. When projecting to a line in 3D though, it's no longer that simple, as $ \mathbf{B} $ will now have two columns. The good news is that we won't have to build that matrix $ \mathbf{B} $ at all - we'll get it cancelled instead. But for now, let's assume that we have it.

Now, we can represent any vector $ \mathbf{b} $ as a linear combination of columns of $ \mathbf{A} $ and $ \mathbf{B} $:

$$
\mathbf{Ax} + \mathbf{By} = \mathbf{b} \tag{1}
$$

We could solve the above equation by introducing a matrix $ \mathbf{M} $ and a vector $ \mathbf{z} $, where:

$$
\begin{align}
\mathbf{M} &= \left[
\begin{array}{c|c}
\mathbf{A} & \mathbf{B}
\end{array}
\right] \notag \\
\mathbf{z} &=
\begin{bmatrix}
\mathbf{x} \\ \mathbf{y}
\end{bmatrix} \notag \\
\end{align}
$$

Then, solving $ \mathbf{Mz} = \mathbf{b} $ for $ \mathbf{z} $ gives us $ \mathbf{x} $ and $ \mathbf{y} $, the latter of which we discard.

However, we are going to do better than that - we are going to left-multiply both sides of (1) by something that will cancel out the $ \mathbf{By} $ term. That something happens to be $ \mathbf{A}^\top $.

That gives us the following equation:

$$
\mathbf{A}^\top \mathbf{Ax} + \mathbf{A}^\top \mathbf{By} = \mathbf{A}^\top \mathbf{b}
$$

I claim that $ \mathbf{A}^\top \mathbf{B} $ is a matrix of zeros and thus $ \mathbf{A^\top By} $ is a vector of zeros for any $ \mathbf{y} $.

What are the elements of $ \mathbf{A}^\top \mathbf{B} $? They are dot products between a column of $ \mathbf{B} $ and a row of $ \mathbf{A}^\top $ (that is a column of $ \mathbf{A} $). Recall that the column spaces of $ \mathbf{A} $ and $ \mathbf{B} $ are orthogonal - so, every column of $ \mathbf{A} $ is orthogonal to every column of $ \mathbf{B} $. Dot products of orthogonal vectors are zeros, so $ \mathbf{A}^\top \mathbf{B} $ is indeed a matrix of zeros. Eliminating $ \mathbf{A}^\top \mathbf{By} $ gives us the well known *Normal Equation*:

$$
\mathbf{A}^\top \mathbf{Ax} = \mathbf{A}^\top \mathbf{b}
$$

Here, $ \mathbf{A}^\top \mathbf{A} $ is a square matrix of the same rank as $ \mathbf{A} $, so this system can easily be solved.

The solution $ \mathbf{x} $ gives us the linear coefficients for the vector in the column space of $ \mathbf{A} $ that's closest (in the Euclidean sense) to the vector $ \mathbf{b} $.

This completes the derivation.

## Finding the Closest Points on Two Lines

Now, let's apply our new knowledge to solve a different problem:

> You have two lines in 3D. Let's call them $ \mathbf{L_1} $ and $ \mathbf{L_2} $. Each line is defined by a point on the line and a vector of arbitrary (nonzero) length, defining the line's direction. So, we've got points $ \mathbf{p_1} $ and $ \mathbf{p_2} $ and vectors $ \mathbf{v_1} $ and $ \mathbf{v_2} $. Your task is to find a point on each line (let's call them $ \mathbf{k_1} $ and $ \mathbf{k_2}) $ such that the distance between those two points is minimised.

Let's try to solve that problem with Linear Algebra. The key to the solution is to realize that the line connecting $ \mathbf{k_1} $ and $ \mathbf{k_2} $ is going to be orthogonal to both $ \mathbf{L_1} $ and $ \mathbf{L_2} $. Let's give the vector from $ \mathbf{k_1} $ to $ \mathbf{k_2} $ a name. Let's call it $ \mathbf{u} $. Now, we can write the following equation:

$$
\mathbf{p_1} + x_1 \mathbf{v_1} + \mathbf{u} = \mathbf{p_2} + x_2 \mathbf{v_2}
$$

There, $ x_1 $ and $ x_2 $ are arbitrary scalars that move us from $ \mathbf{p_i} $ along $ \mathbf{v_i} $.

Now, let's introduce some auxiliary values:

$$
\begin{align}
\mathbf{A} &= \left[
\begin{array}{c|c}
\mathbf{v_1} & -\mathbf{v_2}
\end{array}
\right] \notag \\
\mathbf{b} &= \mathbf{p_2} - \mathbf{p_1} \notag \\
\mathbf{x} &=
\begin{bmatrix}
x_1 \\ x_2
\end{bmatrix} \notag \\
\end{align}
$$

Those values allow us to rewrite the equation above in a matrix form:
$$
\mathbf{Ax} + \mathbf{u} = \mathbf{b}
$$

Now, let's use the same trick we did before to get rid of $ \mathbf{u} $:

$$
\mathbf{A}^\top \mathbf{Ax} + \mathbf{A}^\top \mathbf{u} = \mathbf{A}^\top \mathbf{b}
$$

$ \mathbf{A}^\top \mathbf{u} $ is a vector of zeros because $ \mathbf{u} $ is orthogonal to all columns of $ \mathbf{A} $. That elimination produces the normal equation again:

$$
\mathbf{A}^\top \mathbf{Ax} = \mathbf{A}^\top \mathbf{b}
$$

So, this problem can also be solved by a method of least squares. Let's just do it:

```c++
std::pair<glm::vec3, glm::vec3> closestPointsOnTwoLines(
    glm::vec3 p1, glm::vec3 v1, glm::vec3 p2, glm::vec3 v2)
{
    // Note that in glm, a 2x3 matrix means a 2-column, 3-row, matrix.
    // In Linear Algebra, we would have called such a matrix a 3x2 one.
    const glm::mat2x3 A(v1, -v2);

    const auto b = p2 - p1;
    const auto At(glm::transpose(A));
    const auto invAtA(glm::inverse(At * A));
    const auto x = invAtA * At * b;

    return { p1 + v1 * x[0], p2 + v2 * x[1] };
}
```

## Closing Words

The method of least squares seems to be little known among the software developers and yet is quite useful in both 2D and 3D graphics and beyond. This article should have given you all the knowledge you need to apply it to your projects.

BTW, the [author](/about/) may be available for contracting work.