package main

import "fmt"
import "math"

type vector2 struct {
    x float64
    y float64
}

type particle struct {
    position         vector2
    velocity         vector2
    previousPosition vector2
}

type constants struct {
    g        vector2

    epsilon  float64
    damping  float64
    h        float64
    k        float64
    kNear    float64
    density0 float64
    timeStep float64
    sigma    float64
    beta     float64
}

type boundingBox struct {
    top    float64
    left   float64
    bottom float64
    right  float64
}

func (bb * boundingBox) collide(p * particle, c constants) {

    if p.position.y < bb.bottom {

        p.position.y = bb.bottom + c.epsilon
        p.velocity.y = -p.velocity.y * (1 - c.damping)

    } else if p.position.y > bb.top {

        p.position.y = bb.top - c.epsilon
        p.velocity.y = -p.velocity.y * (1 - c.damping)

    } else if p.position.x < bb.left {

        p.position.x = bb.left + c.epsilon
        p.velocity.y = -p.velocity.x * (1 - c.damping)

    } else if p.position.x > bb.right {

        p.position.x = bb.right - c.epsilon
        p.velocity.y = -p.velocity.x * (1 - c.damping)

    }
}

func (v * vector2) mag() float64 {
    return math.Sqrt(v.x * v.x + v.y * v.y)
}

func (v1 * vector2) multiply(v2 vector2) vector2 {
    return vector2{x : v1.x * v2.x, y : v1.y * v2.y}
}

func (v1 * vector2) divide(v2 vector2) vector2 {
    return vector2{x : v1.x / v2.x, y : v1.y / v2.y}
}

func (v1 * vector2) multiplyF(p float64) vector2 {
    return vector2{x : v1.x * p, y : v1.y * p}
}

func (v1 * vector2) divideF(p float64) vector2 {
    return vector2{x : v1.x / p, y : v1.y / p}
}

func (v1 * vector2) minus(v2 vector2) vector2 {
    return vector2{x : v1.x - v2.x, y : v1.y - v2.y}
}

func (v1 * vector2) plus(v2 vector2) vector2 {
    return vector2{x : v1.x + v2.x, y : v1.y + v2.y}
}

func (v1 * vector2) subtract(v2 vector2) {
    v1.x -= v2.x
    v1.y -= v2.y
}

func (v1 * vector2) add(v2 vector2) {
    v1.x += v2.x
    v1.y += v2.y
}

func (v1 * vector2) subtractF(s float64) {
    v1.x -= s
    v1.y -= s
}

func (v1 * vector2) addF(s float64) {
    v1.x += s
    v1.y += s
}


func distSq(p vector2, p2 vector2) float64 {
    p.subtract(p2)
    return math.Pow(p.x, 2) + math.Pow(p.y, 2)
}

func dist(p1 vector2, p2 vector2) float64 {
    return math.Sqrt( distSq(p1, p2) )
}

func unitVec(p1 vector2, p2 vector2) vector2 {
    d := p1.minus(p2)
    return d.divideF(dist(p1, p2))
}

func (p * particle) updateVelocity(timeStep float64) {
    diffPos := p.position.minus(p.previousPosition)
    p.velocity = diffPos.divideF(timeStep)
}

func (p * particle) integrate(timeStep float64) {
    p.previousPosition = p.position
    p.position.add(p.velocity.multiplyF(timeStep))
}

func applyExternalForces(p * particle, c constants) {

    // Gravity
    p.velocity.add(c.g.multiplyF(c.timeStep))
}

func main() {

    point1 := vector2{x : 3.0, y : 2.0}
    point2 := vector2{x : 12.0, y : 5.0}

    fmt.Println(point1.mag())
    fmt.Println(dist(point1, point2))
    fmt.Println(point1.minus(point2))
    fmt.Println(unitVec(point1, point2))
}
