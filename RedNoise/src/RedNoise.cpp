#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <CanvasPoint.h>
#include <Colour.h>
#include <CanvasTriangle.h>
#include <ModelTriangle.h>
#include <glm/glm.hpp>
#include <RayTriangleIntersection.h>
#include <TexturePoint.h>
#include <TextureMap.h>
#include <array>
#include <chrono>
#include <thread>

#define WIDTH 640
#define HEIGHT 480

CanvasPoint getCanvasIntersectionPoint(glm::vec3 cameraPosition, glm::mat3 cameraOrientation, glm::vec3 vertexPosition, float focalLength, float imageScale) {
	float x = cameraPosition.x - vertexPosition.x;
	float y = cameraPosition.y - vertexPosition.y;
	float z = cameraPosition.z - vertexPosition.z;
	glm::vec3 cameraToPoints = glm::vec3(x, y, z);
	glm::vec3 corrected = cameraToPoints * cameraOrientation;
	float u = (focalLength * (corrected.x / corrected.z)) * (-imageScale) + (WIDTH / 2);
	float v = (focalLength * (corrected.y / corrected.z)) * imageScale + (HEIGHT / 2);
	CanvasPoint cPoint = CanvasPoint::CanvasPoint(std::round(u), std::round(v));
	cPoint.depth = corrected.z;
	return cPoint;
}


glm::vec3 get3DIntersectionPoint(glm::vec3 cameraPosition, glm::mat3 cameraOrientation, CanvasPoint cPoint, float focalLength, float imageScale) {
	float u = cPoint.x;
	float v = cPoint.y;
	float xOverz = (u - (WIDTH / 2)) / (focalLength * (imageScale));
	float yOverz = (v - (HEIGHT / 2)) / (focalLength * (-imageScale));
	float xOvery = xOverz / yOverz;
	float z = cameraPosition.z - focalLength;
	float x = xOverz * z;
	float y = yOverz * z;
	glm::vec3 corrected = glm::vec3(x, y, z);
	glm::vec3 vertexPosition = (cameraOrientation) * corrected;
	return vertexPosition;
}

RayTriangleIntersection getClosestIntersection(glm::vec3 cameraPosition, glm::vec3 rayDirection, std::vector<ModelTriangle>& models) {
	RayTriangleIntersection closest = RayTriangleIntersection(glm::vec3(), 999999, ModelTriangle(), 0);
	glm::vec3 solution;
	for (int i = 0; i < models.size(); i++) {
		glm::vec3 e0 = models[i].vertices[1] - models[i].vertices[0];
		glm::vec3 e1 = models[i].vertices[2] - models[i].vertices[0];
		glm::vec3 SPVector = cameraPosition - models[i].vertices[0];
		glm::mat3 DEMatrix(-rayDirection, e0, e1);
		glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;
		float u = possibleSolution.y;
		float v = possibleSolution.z;
		if (possibleSolution.x < closest.distanceFromCamera && (u >= 0.0) && (u <= 1.0) && (v >= 0.0) && (v <= 1.0) && (u + v) <= 1.0 && possibleSolution.x > 0) {
				glm::vec3 point = cameraPosition + rayDirection * possibleSolution.x;
				closest = RayTriangleIntersection(point, possibleSolution.x, models[i], i);
				closest.triangleIndex = i;				
		}
	}
	return closest;
}

RayTriangleIntersection getClosestShadow(glm::vec3 cameraPosition, glm::vec3 rayDirection, std::vector<ModelTriangle>& models) {
	RayTriangleIntersection closest = RayTriangleIntersection(glm::vec3(), 999999, ModelTriangle(), 0);
	glm::vec3 solution;
	for (int i = 0; i < models.size(); i++) {
		glm::vec3 e0 = models[i].vertices[1] - models[i].vertices[0];
		glm::vec3 e1 = models[i].vertices[2] - models[i].vertices[0];
		glm::vec3 SPVector = cameraPosition - models[i].vertices[0];
		glm::mat3 DEMatrix(-rayDirection, e0, e1);
		glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;
		float u = possibleSolution.y;
		float v = possibleSolution.z;
		if (possibleSolution.x < closest.distanceFromCamera && (u >= 0.0) && (u <= 1.0) && (v >= 0.0) && (v <= 1.0) && (u + v) <= 1.0 && possibleSolution.x > 0 && models[i].colour.name != "Glass") {
			glm::vec3 point = cameraPosition + rayDirection * possibleSolution.x;
			closest = RayTriangleIntersection(point, possibleSolution.x, models[i], i);
			closest.triangleIndex = i;
		}
	}
	return closest;
}


RayTriangleIntersection getClosestReflection(glm::vec3 cameraPosition, glm::vec3 rayDirection, std::vector<ModelTriangle>& models) {
	RayTriangleIntersection closest = RayTriangleIntersection(glm::vec3(), 999999, ModelTriangle(), 0);
	glm::vec3 solution;
	for (int i = 0; i < models.size(); i++) {
		glm::vec3 e0 = models[i].vertices[1] - models[i].vertices[0];
		glm::vec3 e1 = models[i].vertices[2] - models[i].vertices[0];
		glm::vec3 SPVector = cameraPosition - models[i].vertices[0];
		glm::mat3 DEMatrix(-rayDirection, e0, e1);
		glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;
		float u = possibleSolution.y;
		float v = possibleSolution.z;
		if (possibleSolution.x < closest.distanceFromCamera && (u >= 0.0) && (u <= 1.0) && (v >= 0.0) && (v <= 1.0) && (u + v) <= 1.0 && possibleSolution.x > 0 && models[i].colour.name != "Mirror") {
			glm::vec3 point = cameraPosition + rayDirection * possibleSolution.x;
			closest = RayTriangleIntersection(point, possibleSolution.x, models[i], i);
			closest.triangleIndex = i;
		}
	}
	return closest;
}

RayTriangleIntersection getClosestRefraction(glm::vec3 cameraPosition, glm::vec3 rayDirection, std::vector<ModelTriangle>& models) {
	RayTriangleIntersection closest = RayTriangleIntersection(glm::vec3(), 999999, ModelTriangle(), 0);
	glm::vec3 solution;
	for (int i = 0; i < models.size(); i++) {
		glm::vec3 e0 = models[i].vertices[1] - models[i].vertices[0];
		glm::vec3 e1 = models[i].vertices[2] - models[i].vertices[0];
		glm::vec3 SPVector = cameraPosition - models[i].vertices[0];
		glm::mat3 DEMatrix(-rayDirection, e0, e1);
		glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;
		float u = possibleSolution.y;
		float v = possibleSolution.z;
		if (possibleSolution.x < closest.distanceFromCamera && (u >= 0.0) && (u <= 1.0) && (v >= 0.0) && (v <= 1.0) && (u + v) <= 1.0 && possibleSolution.x > 0 && models[i].colour.name != "Glass") {
			glm::vec3 point = cameraPosition + rayDirection * possibleSolution.x;
			closest = RayTriangleIntersection(point, possibleSolution.x, models[i], i);
			closest.triangleIndex = i;
		}
	}
	return closest;
}

glm::vec3 orbit(glm::vec3 cameraPosition, float theta) {
	glm::vec3 newPosition = cameraPosition;
	glm::mat3 rotation = glm::mat3(cos(theta), 0, sin(theta), 0, 1, 0, -sin(theta), 0, cos(theta));
	newPosition = rotation * cameraPosition;
	return newPosition;
}

glm::vec3 cross(glm::vec3 a, glm::vec3 b) {
	glm::vec3 result = glm::vec3(0, 0, 0);
	result[0] = (a[1] * b[2]) - (a[2] * b[1]);
	result[1] = (a[2] * b[0]) - (a[0] * b[2]);
	result[2] = (a[0] * b[1]) - (a[1] * b[0]);
	return result;
}

glm::mat3 lookAt(glm::vec3 center, glm::vec3 cameraPosition, glm::mat3 orientation) {
	glm::vec3 forward = glm::normalize((cameraPosition - center));
	glm::vec3 right = cross(glm::vec3(0, 1, 0), forward);
	glm::vec3 up = cross(forward, right);
	glm::mat3 newOrientation = glm::mat3(right.x, right.y, right.z, up.x, up.y, up.z, forward.x, forward.y, forward.z);
	return newOrientation;
}

glm::vec3 moveCamera(glm::vec3 cameraPosition, glm::mat3 cameraOrientation, SDL_Event event) {
	float offset = .01;
	float theta = .02;
	glm::vec3 newPosition = cameraPosition;
	if (event.key.keysym.sym == SDLK_LEFT) {
		newPosition.x += offset;
	}
	else if (event.key.keysym.sym == SDLK_RIGHT) {
		newPosition.x -= offset;
	}
	else if (event.key.keysym.sym == SDLK_UP) {
		newPosition.y -= offset;
	}
	else if (event.key.keysym.sym == SDLK_DOWN) {
		newPosition.y += offset;
	}
	else if (event.key.keysym.sym == SDLK_f) {
		newPosition.z -= offset;

	}
	else if (event.key.keysym.sym == SDLK_b) {
		newPosition.z += offset;
	}
	else if (event.key.keysym.sym == SDLK_a) {
		glm::mat3 rotation = glm::mat3(cos(theta), 0, -sin(theta), 0, 1, 0, sin(theta), 0, cos(theta));
		newPosition = rotation * cameraPosition;
	}
	else if (event.key.keysym.sym == SDLK_d) {
		glm::mat3 rotation = glm::mat3(cos(theta), 0, sin(theta), 0, 1, 0, -sin(theta), 0, cos(theta));
		newPosition = rotation * cameraPosition;
	}
	else if (event.key.keysym.sym == SDLK_s) {
		glm::mat3 rotation = glm::mat3(1, 0, 0, 0, cos(theta), -sin(theta), 0, sin(theta), cos(theta));
		newPosition = rotation * cameraPosition;
	}
	else if (event.key.keysym.sym == SDLK_w) {
		glm::mat3 rotation = glm::mat3(1, 0, 0, 0, cos(theta), sin(theta), 0, -sin(theta), cos(theta));
		newPosition = rotation * cameraPosition;
	}
	return newPosition;
}

glm::vec3 moveLight(glm::vec3 lightSource, SDL_Event event) {
	float offset = .1;
	glm::vec3 newPosition = lightSource;
	if (event.key.keysym.sym == SDLK_j) {
		newPosition.x -= offset;
	}
	else if (event.key.keysym.sym == SDLK_l) {
		newPosition.x += offset;
	}
	else if (event.key.keysym.sym == SDLK_i) {
		newPosition.y += offset;
	}
	else if (event.key.keysym.sym == SDLK_k) {
		newPosition.y -= offset;
	}
	else if (event.key.keysym.sym == SDLK_u) {
		newPosition.z -= offset;

	}
	else if (event.key.keysym.sym == SDLK_o) {
		newPosition.z += offset;
	}
	return newPosition;
}


glm::mat3 rotateCamera(glm::mat3 cameraOrientation, SDL_Event event) {
	float theta = .02;
	glm::mat3 newOrientation = cameraOrientation;
	if (event.key.keysym.sym == SDLK_a) {
	glm::mat3 rotation = glm::mat3(cos(theta), 0, -sin(theta), 0, 1, 0, sin(theta), 0, cos(theta));
	newOrientation = rotation * cameraOrientation;
	}
	else if (event.key.keysym.sym == SDLK_d) {
	glm::mat3 rotation = glm::mat3(cos(theta), 0, sin(theta), 0, 1, 0, -sin(theta), 0, cos(theta));
	newOrientation = rotation * cameraOrientation;
	}
	else if (event.key.keysym.sym == SDLK_s) {
	glm::mat3 rotation = glm::mat3(1, 0, 0, 0, cos(theta), -sin(theta), 0, sin(theta), cos(theta));
	newOrientation = rotation * cameraOrientation;
	}
	else if (event.key.keysym.sym == SDLK_w) {
	glm::mat3 rotation = glm::mat3(1, 0, 0, 0, cos(theta), sin(theta), 0, -sin(theta), cos(theta));
	newOrientation = rotation * cameraOrientation;
	}
	return newOrientation;
}

float containAngle(float angle) {
	if (angle < 0) return 0;
	else if (angle > 1) return 1.0;
	else return angle;
}

CanvasTriangle sort(CanvasTriangle triangle) {
	CanvasTriangle sorted;
	if (triangle.v0().y <= triangle.v1().y && triangle.v0().y <= triangle.v2().y && triangle.v1().y < triangle.v2().y) {
		sorted = triangle;
	}
	else if (triangle.v0().y <= triangle.v1().y && triangle.v0().y <= triangle.v2().y && triangle.v1().y >= triangle.v2().y) {
		sorted.v0() = triangle.v0();
		sorted.v1() = triangle.v2();
		sorted.v2() = triangle.v1();
	}
	else if (triangle.v1().y <= triangle.v2().y && triangle.v1().y <= triangle.v0().y && triangle.v2().y < triangle.v0().y) {
		sorted.v0() = triangle.v1();
		sorted.v1() = triangle.v2();
		sorted.v2() = triangle.v0();
	}
	else if (triangle.v1().y <= triangle.v2().y && triangle.v1().y <= triangle.v0().y && triangle.v2().y >= triangle.v0().y) {
		sorted.v0() = triangle.v1();
		sorted.v1() = triangle.v0();
		sorted.v2() = triangle.v2();
	}
	else if (triangle.v2().y <= triangle.v1().y && triangle.v2().y <= triangle.v0().y && triangle.v1().y < triangle.v0().y) {
		sorted.v0() = triangle.v2();
		sorted.v1() = triangle.v1();
		sorted.v2() = triangle.v0();
	}
	else if (triangle.v2().y <= triangle.v1().y && triangle.v2().y <= triangle.v0().y && triangle.v1().y >= triangle.v0().y) {
		sorted.v0() = triangle.v2();
		sorted.v1() = triangle.v0();
		sorted.v2() = triangle.v1();
	}
	return sorted;
}

void draw(DrawingWindow &window) {
	window.clearPixels();
	for (size_t y = 0; y < window.height; y++) {
		for (size_t x = 0; x < window.width; x++) {
			float red = rand() % 256;
			float green = 0.0;
			float blue = 0.0;
			uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
			window.setPixelColour(x, y, colour);
		}
	}
}

glm::vec3 vNormalCalculator(ModelTriangle triangle, glm::vec3 vertex, std::vector<ModelTriangle> models) {
	glm::vec3 vNormal = triangle.normal;
	for (int i = 0; i < models.size(); i++) {
		if (models[i].vertices[0] == vertex || models[i].vertices[1] == vertex || models[i].vertices[2] == vertex) {
			vNormal += models[i].normal;
		}
	}
	return glm::normalize(vNormal);
}

std::vector<glm::vec3> multiLight(glm::vec3 center, int radius, float spread) {
	float x = center.x;
	float y = center.y;
	float z = center.z;
	std::vector<glm::vec3> lightPoints;
	//lightPoints.push_back(center);
	for (int i = -radius; i <= radius; i++) {
		for (int j = -radius; j <= radius; j++) {
			if (std::fabs(j) + std::fabs(i) <= radius)  lightPoints.push_back(glm::vec3(x + i * spread, y, z + j * spread));
		}
	}
	return lightPoints;
}

void rayTrace(DrawingWindow& window, glm::vec3 cameraPosition, glm::mat3 cameraOrientation, glm::vec3 lightSource, float focalLength, float imageScale, int shadingCode, int shadowCode, std::vector<ModelTriangle> models) {
	window.clearPixels();
	std::vector<glm::vec3> lightPoints;
	if (shadowCode == 2) lightPoints = multiLight(lightSource, 3, .05);
	for (size_t y = 0; y < window.height; y++) {
		for (size_t x = 0; x < window.width; x++) {
			CanvasPoint cPoint = CanvasPoint::CanvasPoint(x, y, focalLength);
			glm::vec3 threePoint = get3DIntersectionPoint(cameraPosition, cameraOrientation, cPoint, focalLength, imageScale);
			glm::vec3 rayDirection = threePoint - cameraPosition;
			RayTriangleIntersection  intersection = getClosestIntersection(cameraPosition, rayDirection, models);
			glm::vec3 lightRay = glm::normalize(intersection.intersectionPoint - lightSource);
			glm::vec3 pointToLight = glm::normalize(lightSource - intersection.intersectionPoint);
			glm::vec3 view = glm::normalize(intersection.intersectionPoint - cameraPosition);

			//Lighting Settings
			float pL = 0;
			float pLPower = 1;
			float sLGloss = 256;
			double sLScale = 1.0;

			//Vertext Normals
			if (shadingCode == 2 || shadingCode == 3) {
				glm::vec3 A = intersection.intersectedTriangle.vertices[0];
				glm::vec3 B = intersection.intersectedTriangle.vertices[1];
				glm::vec3 C = intersection.intersectedTriangle.vertices[2];
				glm::vec3 P = intersection.intersectionPoint;

				glm::vec3 n0 = vNormalCalculator(intersection.intersectedTriangle, A, models);
				glm::vec3 n1 = vNormalCalculator(intersection.intersectedTriangle, B, models);
				glm::vec3 n2 = vNormalCalculator(intersection.intersectedTriangle, C, models);

				//float ABCArea = glm::length(glm::cross((A - B), (A - C)));
				/*float u = glm::length(glm::cross((C - A), (C - P))) / ABCArea;
				float v = glm::length(glm::cross((A - B), (A - P))) / ABCArea;
				float w = glm::length(glm::cross((B - C), (B - P))) / ABCArea;*/

				float det = ((B.y - C.y) * (A.x - C.x)) + ((C.x - B.x) * (A.y - C.y));

				float u = (((B.y - C.y) * (P.x - C.x)) + ((C.x - B.x) * (P.y - C.y))) / det;
				float v = (((C.y - A.y) * (P.x - C.x)) + ((A.x - C.x) * (P.y - C.y))) / det;
				float w = 1 - u - v;

				glm::vec3 pointNormal = u * n0 + v * n1 + w * n2;

				//int shadingCode = 3; // 1 is flat shading, 2 is gouraud, 3 is phong

				if (shadingCode == 2) {
					//Incident Lightning
					float incidentAngle1 = glm::dot(n0, pointToLight);
					float incidentAngle2 = glm::dot(n1, pointToLight);
					float incidentAngle3 = glm::dot(n2, pointToLight);
					//Spherical Lighting
					glm::vec3 reflectionVector1 = glm::normalize(lightRay) - 2.0f * n0 * glm::dot(glm::normalize(lightRay), n0);
					glm::vec3 reflectionVector2 = glm::normalize(lightRay) - 2.0f * n1 * glm::dot(glm::normalize(lightRay), n1);
					glm::vec3 reflectionVector3 = glm::normalize(lightRay) - 2.0f * n2 * glm::dot(glm::normalize(lightRay), n2);
					float sL1 = std::fabs(std::pow(glm::dot(view, reflectionVector1), sLGloss) * sLScale);
					float sL2 = std::fabs(std::pow(glm::dot(view, reflectionVector2), sLGloss) * sLScale);
					float sL3 = std::fabs(std::pow(glm::dot(view, reflectionVector3), sLGloss) * sLScale);
					//Proximity Lighting
					float pL1 = ((pLPower) / (4 * 3.1415 * glm::length(lightRay) * glm::length(lightRay))) + incidentAngle1 + sL1;
					float pL2 = ((pLPower) / (4 * 3.1415 * glm::length(lightRay) * glm::length(lightRay))) + incidentAngle2 + sL2;
					float pL3 = ((pLPower) / (4 * 3.1415 * glm::length(lightRay) * glm::length(lightRay))) + incidentAngle3 + sL3;
					pL = u * pL1 + v * pL2 + w * pL3;
				}
				
				//Phong Shading
				else if (shadingCode == 3) {
					float incidentAngle = glm::dot(pointNormal, pointToLight);
					glm::vec3 reflectionVectorP = glm::normalize(lightRay) - 2.0f * pointNormal * glm::dot(glm::normalize(lightRay), pointNormal);
					float sL = std::fabs(std::pow(glm::dot(view, reflectionVectorP), sLGloss) * sLScale);
					pL = ((pLPower) / (4 * 3.1415 * glm::length(lightRay) * glm::length(lightRay))) + incidentAngle + sL;
				}
			}
			//Flat Shading
			if (pL != pL || shadingCode == 1) {
				float incidentAngle = glm::dot(intersection.intersectedTriangle.normal, pointToLight);
				glm::vec3 reflectionVectorF = glm::normalize(lightRay) - 2.0f * glm::normalize(intersection.intersectedTriangle.normal) * glm::dot(glm::normalize(lightRay), glm::normalize(intersection.intersectedTriangle.normal));
				float sL = std::fabs(std::pow(glm::dot(view, reflectionVectorF), sLGloss) * sLScale);
				pL = ((pLPower) / (4 * 3.1415 * glm::length(lightRay) * glm::length(lightRay))) + incidentAngle + sL;
			}
			//Ambient Lightning
			float ambient = .6;
			if (pL < ambient) pL = ambient;
			
			//Reflection
			glm::vec3 normal = glm::normalize(intersection.intersectedTriangle.normal);
			//glm::vec3 reflectionVector = glm::normalize(lightRay) - 2.0f * pointNormal * glm::dot(glm::normalize(lightRay), pointNormal);
			//glm::vec3 reflectionVector = glm::normalize(lightRay) - 2.0f * glm::normalize(intersection.intersectedTriangle.normal) * glm::dot(glm::normalize(lightRay), glm::normalize(intersection.intersectedTriangle.normal));
			glm::vec3 reflectionVector = glm::normalize(view - (2.0f * (glm::dot(view, normal)) * normal));
			Colour col = intersection.intersectedTriangle.colour;
			float mirrorMod = 1.0;
			if (intersection.intersectedTriangle.colour.name == "Mirror") {
				//glm::vec3 bias = normal * 1.0f;
				RayTriangleIntersection reflection = getClosestReflection(intersection.intersectionPoint, reflectionVector, models);
				col = reflection.intersectedTriangle.colour;
				//std::cout << col << std::endl;
				mirrorMod = .7;
			}

			//Refraction 
			if (intersection.intersectedTriangle.colour.name == "Glass") {
				float glassRI = 1.5;
				float airRI = 1.0;
				float N1 = airRI;
				float N2 = glassRI;
				float incidentAngle = glm::dot(normal, view);
				float RI = N1 / N2;
				float theta1 = incidentAngle;
				glm::vec3 I = view;
				float c1 = glm::dot(normal, I);
				float refractionAngle = std::asin(std::sin(theta1) * RI);
				float theta2 = refractionAngle;
				float c2 = (float)std::sqrt(1 - (std::pow(RI, 2) * (1 - std::pow(std::cos(theta1), 2))));
				glm::vec3 tRay = RI * I + (RI * c1 - c2) * normal;
				glm::vec3 bias = normal * 0.0f;
				RayTriangleIntersection refraction = getClosestRefraction(intersection.intersectionPoint - bias, tRay, models);
				RayTriangleIntersection reflection = getClosestRefraction(intersection.intersectionPoint + bias, reflectionVector, models);
				//float refractionAngle = glm::dot(tRay, normal);
				float F1 = std::pow((N2 * std::cos(theta1) - N1 * std::cos(theta2)) / (N2 * std::cos(theta1) + N1 * std::cos(theta2)), 2);
				float F2 = std::pow((N1 * std::cos(theta2) - N2 * std::cos(theta1)) / (N1 * std::cos(theta2) + N2 * std::cos(theta1)), 2);
				float FReflect = (F1 * F2) / 2;
				//std::cout << FReflect << std::endl;
				float FRefract = 1 - FReflect;
				Colour reflectC = reflection.intersectedTriangle.colour;
				Colour refractC = refraction.intersectedTriangle.colour;
				col.red = FReflect * reflectC.red + FRefract * refractC.red;
				col.green = FReflect * reflectC.green+ FRefract * refractC.green;
				col.blue = FReflect * reflectC.blue + FRefract * refractC.blue;
			}

			//Colouring
			float red = std::min(col.red * pL * mirrorMod, 255.0f);
			float green = std::min(col.green * pL * mirrorMod, 255.0f);
			float blue = std::min(col.blue * pL * mirrorMod, 255.0f);
			//Shadow Calculation
			//Hard Shadows
			if (shadowCode == 1) {
				RayTriangleIntersection  shadowIntersection = getClosestShadow(lightSource, lightRay, models);
				float shadowValue = ambient;
				if ((int)intersection.triangleIndex != (int)shadowIntersection.triangleIndex) {
					red = red * shadowValue;
					green = green * shadowValue;
					blue = blue * shadowValue;
				}
			}
			//Soft Shadows
			if (shadowCode == 2) {
				RayTriangleIntersection  shadowIntersection = getClosestShadow(lightSource, lightRay, models);
				float shadowValue = 1;
				//if ((int)intersection.triangleIndex != (int)shadowIntersection.triangleIndex) {
					for (int i = 0; i < lightPoints.size(); i++) {
						glm::vec3 lRay = intersection.intersectionPoint - lightPoints[i];
						RayTriangleIntersection  multiShadowIntersection = getClosestShadow(lightPoints[i], lRay, models);
						if ((int)intersection.triangleIndex != (int)multiShadowIntersection.triangleIndex) shadowValue -= .02;
					}
					shadowValue = std::max(ambient, shadowValue);
					red = red * shadowValue;
					green = green * shadowValue;
					blue = blue * shadowValue;
				//}
			}
			uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
			window.setPixelColour(x, y, colour);
		}
	}
	/*glm::vec3 v = lightSource - cameraPosition;
	v = v * cameraOrientation;
	float ui = (focalLength * (-v.x * imageScale / v.z) + (WIDTH / 2));
	float vi = (focalLength * (v.y * imageScale / v.z) + (HEIGHT / 2));
	Colour col = Colour::Colour(255, 255, 255);
	float red = col.red;
	float green = col.green;
	float blue = col.blue;
	uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
	window.setPixelColour(ui, vi, colour);
	for (int i = 0; i < lightPoints.size(); i++) {
		glm::vec3 v = lightPoints[i] - cameraPosition;
		v = v * cameraOrientation;
		float ui2 = (focalLength * (-v.x * imageScale / v.z) + (WIDTH / 2));
		float vi2 = (focalLength * (v.y * imageScale / v.z) + (HEIGHT / 2));
		window.setPixelColour(ui2, vi2, colour);
	}*/
}

void rayTraceP(DrawingWindow& window, glm::vec3 cameraPosition, glm::mat3 cameraOrientation, glm::vec3 lightSource, float focalLength, float imageScale, int shadingCode, int shadowCode, std::vector<ModelTriangle> models, int scale, std::vector<glm::vec3> dimensions, std::vector<std::vector<std::vector<float>>> pMap) {
	window.clearPixels();
	float xmin = dimensions[1][0];
	float ymin = dimensions[1][1];
	float zmin = dimensions[1][2];
	std::vector<glm::vec3> lightPoints;
	if (shadowCode == 2) lightPoints = multiLight(lightSource, 3, .05);
	for (size_t y = 0; y < window.height; y++) {
		for (size_t x = 0; x < window.width; x++) {
			CanvasPoint cPoint = CanvasPoint::CanvasPoint(x, y, focalLength);
			glm::vec3 threePoint = get3DIntersectionPoint(cameraPosition, cameraOrientation, cPoint, focalLength, imageScale);
			glm::vec3 rayDirection = threePoint - cameraPosition;
			RayTriangleIntersection  intersection = getClosestIntersection(cameraPosition, rayDirection, models);
			glm::vec3 lightRay = glm::normalize(intersection.intersectionPoint - lightSource);
			glm::vec3 pointToLight = glm::normalize(lightSource - intersection.intersectionPoint);
			glm::vec3 view = glm::normalize(intersection.intersectionPoint - cameraPosition);

			//Lighting Settings
			float pL = 0;
			float pLPower = 1;
			float sLGloss = 256;
			double sLScale = 1.0;
			float minLight = .35;

			//Vertext Normals
			if (shadingCode == 2 || shadingCode == 3) {
				glm::vec3 A = intersection.intersectedTriangle.vertices[0];
				glm::vec3 B = intersection.intersectedTriangle.vertices[1];
				glm::vec3 C = intersection.intersectedTriangle.vertices[2];
				glm::vec3 P = intersection.intersectionPoint;

				glm::vec3 n0 = vNormalCalculator(intersection.intersectedTriangle, A, models);
				glm::vec3 n1 = vNormalCalculator(intersection.intersectedTriangle, B, models);
				glm::vec3 n2 = vNormalCalculator(intersection.intersectedTriangle, C, models);

				float ABCArea = glm::length(glm::cross((A - B), (A - C)));
				/*float u = glm::length(glm::cross((C - A), (C - P))) / ABCArea;
				float v = glm::length(glm::cross((A - B), (A - P))) / ABCArea;
				float w = glm::length(glm::cross((B - C), (B - P))) / ABCArea;*/

				float det = ((B.y - C.y) * (A.x - C.x)) + ((C.x - B.x) * (A.y - C.y));

				float u = (((B.y - C.y) * (P.x - C.x)) + ((C.x - B.x) * (P.y - C.y))) / det;
				float v = (((C.y - A.y) * (P.x - C.x)) + ((A.x - C.x) * (P.y - C.y))) / det;
				float w = 1 - u - v;

				glm::vec3 pointNormal = u * n0 + v * n1 + w * n2;

				//int shadingCode = 3; // 1 is flat shading, 2 is gouraud, 3 is phong

				if (shadingCode == 2) {
					//Incident Lightning
					float incidentAngle1 = glm::dot(n0, pointToLight);
					float incidentAngle2 = glm::dot(n1, pointToLight);
					float incidentAngle3 = glm::dot(n2, pointToLight);
					//Spherical Lighting
					glm::vec3 reflectionVector1 = glm::normalize(lightRay) - 2.0f * n0 * glm::dot(glm::normalize(lightRay), n0);
					glm::vec3 reflectionVector2 = glm::normalize(lightRay) - 2.0f * n1 * glm::dot(glm::normalize(lightRay), n1);
					glm::vec3 reflectionVector3 = glm::normalize(lightRay) - 2.0f * n2 * glm::dot(glm::normalize(lightRay), n2);
					float sL1 = std::fabs(std::pow(glm::dot(view, reflectionVector1), sLGloss) * sLScale);
					float sL2 = std::fabs(std::pow(glm::dot(view, reflectionVector2), sLGloss) * sLScale);
					float sL3 = std::fabs(std::pow(glm::dot(view, reflectionVector3), sLGloss) * sLScale);
					//Proximity Lighting
					float pL1 = ((pLPower) / (4 * 3.1415 * glm::length(lightRay) * glm::length(lightRay))) + incidentAngle1 + sL1;
					float pL2 = ((pLPower) / (4 * 3.1415 * glm::length(lightRay) * glm::length(lightRay))) + incidentAngle2 + sL2;
					float pL3 = ((pLPower) / (4 * 3.1415 * glm::length(lightRay) * glm::length(lightRay))) + incidentAngle3 + sL3;
					pL = u * pL1 + v * pL2 + w * pL3;
				}

				//Phong Shading
				if (shadingCode == 3) {
					float incidentAngle = glm::dot(pointNormal, pointToLight);
					glm::vec3 reflectionVectorP = glm::normalize(lightRay) - 2.0f * pointNormal * glm::dot(glm::normalize(lightRay), pointNormal);
					float sL = std::fabs(std::pow(glm::dot(view, reflectionVectorP), sLGloss) * sLScale);
					pL = ((pLPower) / (4 * 3.1415 * glm::length(lightRay) * glm::length(lightRay))) + incidentAngle + sL;
				}
			}
			//Flat Shading
			if (pL != pL || shadingCode == 1) {
				float incidentAngle = glm::dot(intersection.intersectedTriangle.normal, pointToLight);
				glm::vec3 reflectionVectorF = glm::normalize(lightRay) - 2.0f * glm::normalize(intersection.intersectedTriangle.normal) * glm::dot(glm::normalize(lightRay), glm::normalize(intersection.intersectedTriangle.normal));
				float sL = std::fabs(std::pow(glm::dot(view, reflectionVectorF), sLGloss) * sLScale);
				pL = ((pLPower) / (4 * 3.1415 * glm::length(lightRay) * glm::length(lightRay))) + incidentAngle + sL;
			}
			//Ambient Lightning with Photon Map
			int xPos = std::round((intersection.intersectionPoint[0] - xmin) * scale);
			int yPos = std::round((intersection.intersectionPoint[1] - ymin) * scale);
			int zPos = std::round((intersection.intersectionPoint[2] - zmin) * scale);
			float ambient = pMap[xPos][yPos][zPos];
			pL = pL + ambient;
			if (pL < minLight) pL = minLight;

			//Reflection
			glm::vec3 normal = glm::normalize(intersection.intersectedTriangle.normal);
			//glm::vec3 reflectionVector = glm::normalize(lightRay) - 2.0f * pointNormal * glm::dot(glm::normalize(lightRay), pointNormal);
			//glm::vec3 reflectionVector = glm::normalize(lightRay) - 2.0f * glm::normalize(intersection.intersectedTriangle.normal) * glm::dot(glm::normalize(lightRay), glm::normalize(intersection.intersectedTriangle.normal));
			glm::vec3 reflectionVector = glm::normalize(view - (2.0f * (glm::dot(view, normal)) * normal));
			Colour col = intersection.intersectedTriangle.colour;
			float mirrorMod = 1.0;
			if (intersection.intersectedTriangle.colour.name == "Mirror") {
				//glm::vec3 bias = normal * 1.0f;
				RayTriangleIntersection reflection = getClosestReflection(intersection.intersectionPoint, reflectionVector, models);
				col = reflection.intersectedTriangle.colour;
				//std::cout << col << std::endl;
				mirrorMod = .7;
			}

			if (intersection.intersectedTriangle.colour.name == "Glass") {
				float glassRI = 1.5;
				float airRI = 1.0;
				float N1 = airRI;
				float N2 = glassRI;
				float incidentAngle = glm::dot(normal, view);
				float RI = N1 / N2;
				float theta1 = incidentAngle;
				glm::vec3 I = view;
				float c1 = glm::dot(normal, I);
				float refractionAngle = std::asin(std::sin(theta1) * RI);
				float theta2 = refractionAngle;
				float c2 = (float)std::sqrt(1 - (std::pow(RI, 2) * (1 - std::pow(std::cos(theta1), 2))));
				glm::vec3 tRay = RI * I + (RI * c1 - c2) * normal;
				glm::vec3 bias = normal * 0.0f;
				RayTriangleIntersection refraction = getClosestRefraction(intersection.intersectionPoint - bias, tRay, models);
				RayTriangleIntersection reflection = getClosestRefraction(intersection.intersectionPoint + bias, reflectionVector, models);
				//float refractionAngle = glm::dot(tRay, normal);
				float F1 = std::pow((N2 * std::cos(theta1) - N1 * std::cos(theta2)) / (N2 * std::cos(theta1) + N1 * std::cos(theta2)), 2);
				float F2 = std::pow((N1 * std::cos(theta2) - N2 * std::cos(theta1)) / (N1 * std::cos(theta2) + N2 * std::cos(theta1)), 2);
				float FReflect = (F1 * F2) / 2;
				//std::cout << FReflect << std::endl;
				float FRefract = 1 - FReflect;
				Colour reflectC = reflection.intersectedTriangle.colour;
				Colour refractC = refraction.intersectedTriangle.colour;
				col.red = FReflect * reflectC.red + FRefract * refractC.red;
				col.green = FReflect * reflectC.green + FRefract * refractC.green;
				col.blue = FReflect * reflectC.blue + FRefract * refractC.blue;
			}

			//Colouring
			float red = std::min(col.red * pL * mirrorMod, 255.0f);
			float green = std::min(col.green * pL * mirrorMod, 255.0f);
			float blue = std::min(col.blue * pL * mirrorMod, 255.0f);
			//Shadow Calculation
			//Hard Shadows
			if (shadowCode == 1) {
				RayTriangleIntersection  shadowIntersection = getClosestShadow(lightSource, lightRay, models);
				float shadowValue = .6;
				if ((int)intersection.triangleIndex != (int)shadowIntersection.triangleIndex) {
					red = red * shadowValue;
					green = green * shadowValue;
					blue = blue * shadowValue;
				}
			}
			//Soft Shadows
			if (shadowCode == 2) {
				RayTriangleIntersection  shadowIntersection = getClosestShadow(lightSource, lightRay, models);
				float shadowValue = 1;
				//if ((int)intersection.triangleIndex != (int)shadowIntersection.triangleIndex) {
				for (int i = 0; i < lightPoints.size(); i++) {
					glm::vec3 lRay = intersection.intersectionPoint - lightPoints[i];
					RayTriangleIntersection  multiShadowIntersection = getClosestShadow(lightPoints[i], lRay, models);
					if ((int)intersection.triangleIndex != (int)multiShadowIntersection.triangleIndex) shadowValue -= .02;
				}
				shadowValue = std::max(ambient, shadowValue);
				red = red * shadowValue;
				green = green * shadowValue;
				blue = blue * shadowValue;
				//}
			}
			uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
			window.setPixelColour(x, y, colour);
		}
	}
	/*glm::vec3 v = lightSource - cameraPosition;
	v = v * cameraOrientation;
	float ui = (focalLength * (-v.x * imageScale / v.z) + (WIDTH / 2));
	float vi = (focalLength * (v.y * imageScale / v.z) + (HEIGHT / 2));
	Colour col = Colour::Colour(255, 255, 255);
	float red = col.red;
	float green = col.green;
	float blue = col.blue;
	uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
	window.setPixelColour(ui, vi, colour);
	for (int i = 0; i < lightPoints.size(); i++) {
		glm::vec3 v = lightPoints[i] - cameraPosition;
		v = v * cameraOrientation;
		float ui2 = (focalLength * (-v.x * imageScale / v.z) + (WIDTH / 2));
		float vi2 = (focalLength * (v.y * imageScale / v.z) + (HEIGHT / 2));
		window.setPixelColour(ui2, vi2, colour);
	}*/
}

std::vector<float> interpolateSingleFloats(float from, float to, int numberOfValues) {
	std::vector<float> res;
	float inc = (to - from) / (numberOfValues - 1);
	for (int i = 0; i < numberOfValues; i++) {
		res.push_back(from + i * inc);
	}
	return res;
}

std::vector<float> interpolateSingleFloats2(float one, float two, int numberOfValues) {
	std::vector<float> res;
	float to;
	float from;
	if (one <= two) {
		to = one;
		from = two;
	}
	else {
		to = two;
		from = one;
	}
	int inc = (to - from) / (numberOfValues - 1);
	for (int i = 0; i < numberOfValues; i++) {
		res.push_back(std::round(from + i * inc));
	}

	return res;
}

std::vector<float> interpolateSingleFloats3(float one, float two, int numberOfValues) {
	std::vector<float> res;
	float to;
	float from;
	if (one <= two) {
		to = one;
		from = two;
	}
	else {
		to = two;
		from = one;
	}
	int inc = (to - from) / (numberOfValues - 1);
	for (int i = 0; i < numberOfValues; i++) {
		res.push_back(from + i * inc);
	}

	return res;
}

std::vector<TexturePoint> interpolateTexture(TexturePoint from, TexturePoint to, int numberOfValues) {
	std::vector<TexturePoint> res;
		std::vector<float> one = interpolateSingleFloats(from.x, to.x, numberOfValues);
		std::vector<float> two = interpolateSingleFloats(from.y, to.y, numberOfValues);
		for (int i = 0; i < numberOfValues; i++) {
			TexturePoint step = TexturePoint::TexturePoint(std::round(one[i]), std::round(two[i]));
			res.push_back(step);
		}
		return res;
}

std::vector<CanvasPoint> interpolateThreeElementValues(CanvasPoint from, CanvasPoint to, int numberOfValues) {
	std::vector<CanvasPoint> res;
	std::vector<float> one = interpolateSingleFloats(from.x, to.x, numberOfValues);
	std::vector<float> two = interpolateSingleFloats(from.y, to.y, numberOfValues);
	std::vector<float> three = interpolateSingleFloats(from.depth, to.depth, numberOfValues);
	for (int i = 0; i < numberOfValues; i++) {
		CanvasPoint step = CanvasPoint::CanvasPoint(one[i], two[i], three[i]);
		res.push_back(step);
	}
	return res;
}

void drawLine(DrawingWindow& window, CanvasPoint from, CanvasPoint to, Colour colour) {
	//window.clearPixels();
	float xDiff = to.x - from.x;
	float yDiff = to.y - from.y;
	float numberOfSteps = std::max(std::fabs(xDiff), std::fabs(yDiff));

	std::vector<float> xValues = interpolateSingleFloats(from.x, to.x, std::round(numberOfSteps));
	std::vector<float> yValues = interpolateSingleFloats(from.y, to.y, std::round(numberOfSteps));

	float red = colour.red;
	float green = colour.green;
	float blue = colour.blue;
	uint32_t col = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
	for (int i = 0; i < numberOfSteps; i++) {
		window.setPixelColour(std::round(xValues[i]), std::round(yValues[i]), col);
	}
}


std::vector<glm::vec2> drawLine2(DrawingWindow& window, CanvasPoint from, CanvasPoint to) {
	//window.clearPixels();
	float xDiff = to.x - from.x;
	float yDiff = to.y - from.y;
	float numberOfSteps = std::max(std::fabs(xDiff), std::fabs(yDiff));

	std::vector<float> xValues = interpolateSingleFloats(from.x, to.x, std::round(numberOfSteps));
	std::vector<float> yValues = interpolateSingleFloats(from.y, to.y, std::round(numberOfSteps));

	std::vector<glm::vec2> res;
	for (int i = 0; i < numberOfSteps; i++) {
		glm::vec2 step = glm::vec2(xValues[i], yValues[i]);
		res.push_back(step);
	}
	return res;
}

void drawTriangle(DrawingWindow& window, CanvasTriangle triangle, Colour colour) {
	//window.clearPixels();
	drawLine(window, triangle.v0(), triangle.v1(), colour); 
	drawLine(window, triangle.v1(), triangle.v2(), colour);
	drawLine(window, triangle.v2(), triangle.v0(), colour);
}

void drawFilledTop(DrawingWindow& window, CanvasTriangle triangle, Colour colour, std::vector<std::vector<float>> &depthBuffer) {
	std::vector<glm::vec2> line1 = drawLine2(window, triangle.v0(), triangle.v1());
	std::vector<glm::vec2> line2 = drawLine2(window, triangle.v0(), triangle.v2());
	std::vector<glm::vec2> flat = drawLine2(window, triangle.v1(), triangle.v2());
	int triangleHeight = std::fabs(triangle.v2().y - triangle.v0().y);

	std::vector<float> line1correct;
	std::vector<float> line2correct;

	if (line1.size() == 0 || line2.size() == 0) return;

	/*line1correct = interpolateSingleFloats(line1[0][0], line1[line1.size() - 1][0], triangleHeight+1);
	line2correct = interpolateSingleFloats(line2[0][0], line2[line2.size() - 1][0], triangleHeight+1);

	std::vector<float> dLine1 = interpolateSingleFloats(triangle.v0().depth, triangle.v1().depth, triangleHeight+1);
	std::vector<float> dLine2 = interpolateSingleFloats(triangle.v0().depth, triangle.v2().depth, triangleHeight+1);*/

	std::vector<CanvasPoint> line1All = interpolateThreeElementValues(triangle.v0(), triangle.v1(), triangleHeight + 1);
	std::vector<CanvasPoint> line2All = interpolateThreeElementValues(triangle.v0(), triangle.v2(), triangleHeight + 1);

	for (size_t y = 0; y < triangleHeight+1; y++) {
		/*std::vector<float> strip;
		strip = interpolateSingleFloats3(line1correct[y], line2correct[y], std::fabs(std::round((line1correct[y] - line2correct[y]))));
		std::vector<float> dStrip;
		dStrip = interpolateSingleFloats(dLine1[y], dLine2[y], std::fabs(std::round((line1correct[y] - line2correct[y]))));*/
		std::vector<CanvasPoint> allStrip = interpolateThreeElementValues(line1All[y], line2All[y], std::fabs(std::round((line1All[y].x - line2All[y].x)))*2);

		for (size_t x = 0; x < allStrip.size(); x++) {
			float red = colour.red;
			float green = colour.green;
			float blue = colour.blue;
			uint32_t col = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);

			int yValue = std::round(y + triangle.v0().y);
			int xValue = std::round(allStrip[x].x);
			if (xValue >= 0 && xValue < WIDTH && yValue >= 0 && yValue < HEIGHT) {
				if (depthBuffer[yValue][xValue] < (1 / allStrip[x].depth)) {
					depthBuffer[yValue][xValue] = 1 / allStrip[x].depth;
					window.setPixelColour(xValue, yValue, col);
				}
			}
		}
	}
	//drawLine(window, triangle.v1(), triangle.v2(), colour);
}

void drawFilledBottom(DrawingWindow& window, CanvasTriangle triangle, Colour colour, std::vector<std::vector<float>> &depthBuffer) {
	std::vector<glm::vec2> line1 = drawLine2(window, triangle.v0(), triangle.v1());
	std::vector<glm::vec2> line2 = drawLine2(window, triangle.v0(), triangle.v2());
	std::vector<glm::vec2> flat = drawLine2(window, triangle.v1(), triangle.v2());
	int triangleHeight = std::fabs(triangle.v2().y - triangle.v0().y);

	std::vector<float> line1correct;
	std::vector<float> line2correct;

	if (line1.size() == 0 || line2.size() == 0) return;

	/*line1correct = interpolateSingleFloats(line1[0][0], line1[line1.size() - 1][0], triangleHeight);
	line2correct = interpolateSingleFloats(line2[0][0], line2[line2.size() - 1][0], triangleHeight);

	std::vector<float> dLine1 = interpolateSingleFloats(triangle.v0().depth, triangle.v1().depth, triangleHeight);
	std::vector<float> dLine2 = interpolateSingleFloats(triangle.v0().depth, triangle.v2().depth, triangleHeight);*/

	std::vector<CanvasPoint> line1All = interpolateThreeElementValues(triangle.v0(), triangle.v1(), triangleHeight + 1);
	std::vector<CanvasPoint> line2All = interpolateThreeElementValues(triangle.v0(), triangle.v2(), triangleHeight + 1);

	for (size_t y = 0; y < triangleHeight; y++) {
		/*std::vector<float> strip;
		std::vector<float> dStrip;
		strip = interpolateSingleFloats3(line1correct[y], line2correct[y], std::fabs(std::round((line1correct[y] - line2correct[y]))));
		dStrip = interpolateSingleFloats(dLine1[y], dLine2[y], std::fabs(std::round((line1correct[y] - line2correct[y]))));*/
		std::vector<CanvasPoint> allStrip = interpolateThreeElementValues(line1All[y], line2All[y], std::fabs(std::round((line1All[y].x - line2All[y].x)))*2);
		for (size_t x = 0; x < allStrip.size(); x++) {
			float red = colour.red;
			float green = colour.green;
			float blue = colour.blue;
			uint32_t col = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
			int yValue = std::round(triangle.v0().y - y);
			int xValue = std::round(allStrip[x].x);
			if (xValue >= 0 && xValue < WIDTH && yValue >= 0 && yValue < HEIGHT) {
				if (depthBuffer[yValue][xValue] < (1 / allStrip[x].depth)) {
					depthBuffer[yValue][xValue] = 1 / allStrip[x].depth;
					window.setPixelColour(xValue, yValue, col);
				}
			}
		}
	}
}

std::vector<CanvasPoint> furthestAndNearest(CanvasTriangle triangle) {
	std::vector<CanvasPoint> farNear;
	CanvasPoint farthest = triangle.v0();
	if (triangle.v1().depth > farthest.depth) farthest = triangle.v1();
	if (triangle.v2().depth > farthest.depth) farthest = triangle.v2();
	CanvasPoint nearest = triangle.v0();
	if (triangle.v1().depth < nearest.depth) nearest = triangle.v1();
	if (triangle.v2().depth < nearest.depth) nearest = triangle.v2();
	farNear.push_back(farthest);
	farNear.push_back(nearest);
	return farNear;
}

void drawFilledTextureTop(DrawingWindow& window, CanvasTriangle triangle, Colour colour, std::vector<std::vector<float>> &depthBuffer, TextureMap texture) {
	std::vector<glm::vec2> line1 = drawLine2(window, triangle.v0(), triangle.v1());
	std::vector<glm::vec2> line2 = drawLine2(window, triangle.v0(), triangle.v2());
	std::vector<glm::vec2> flat = drawLine2(window, triangle.v1(), triangle.v2());
	int triangleHeight = std::fabs(triangle.v2().y - triangle.v0().y);

	std::vector<float> line1correct;
	std::vector<float> line2correct;

	if (line1.size() == 0 || line2.size() == 0) return;

	/*line1correct = interpolateSingleFloats(line1[0][0], line1[line1.size() - 1][0], triangleHeight + 1);
	line2correct = interpolateSingleFloats(line2[0][0], line2[line2.size() - 1][0], triangleHeight + 1);

	std::vector<float> dLine1 = interpolateSingleFloats(triangle.v0().depth, triangle.v1().depth, triangleHeight + 1);
	std::vector<float> dLine2 = interpolateSingleFloats(triangle.v0().depth, triangle.v2().depth, triangleHeight + 1);*/

	std::vector<TexturePoint> tLine1 = interpolateTexture(triangle.v0().texturePoint, triangle.v1().texturePoint, triangleHeight + 1);
	std::vector<TexturePoint> tLine2 = interpolateTexture(triangle.v0().texturePoint, triangle.v2().texturePoint, triangleHeight + 1);

	std::vector<CanvasPoint> line1All = interpolateThreeElementValues(triangle.v0(), triangle.v1(), triangleHeight + 1);
	std::vector<CanvasPoint> line2All = interpolateThreeElementValues(triangle.v0(), triangle.v2(), triangleHeight + 1);

	//Texture perspective correction
	std::vector<CanvasPoint> farNear = furthestAndNearest(triangle);
	float z0 = farNear[0].depth;
	float z1 = farNear[1].depth;
	float c0 = farNear[0].texturePoint.y;
	float c1 = farNear[1].texturePoint.y;

	for (size_t y = 0; y < triangleHeight + 1; y++) {
		int size = std::fabs(std::round((line1All[y].x - line2All[y].x))) * 2;
		/*std::vector<float> strip;
		strip = interpolateSingleFloats2(line1correct[y], line2correct[y], size);
		std::vector<float> dStrip;
		dStrip = interpolateSingleFloats(dLine1[y], dLine2[y], size);*/
		std::vector<TexturePoint> tStrip;
		tStrip = interpolateTexture(tLine1[y], tLine2[y], size);
		std::vector<CanvasPoint> allStrip = interpolateThreeElementValues(line1All[y], line2All[y], size);
		int q = triangleHeight - y;

		for (size_t x = 0; x < allStrip.size(); x++) {
			int tIndex = texture.width * std::round(tStrip[x].y) + std::round(tStrip[x].x);
			//int c = (((c0 / z0) * (1 - q)) + ((c1 / z1) * q)) / (((1 / z0) * (1 - q)) + ((1 / z1) * q));
			//std::cout << c << std::endl;
			//int tIndex = texture.width * c + std::round(tStrip[x].x);
			if (tIndex >= texture.pixels.size()) tIndex = texture.pixels.size() - 1;
			uint32_t col = texture.pixels[tIndex];

			int yValue = std::round(y + triangle.v0().y);
			int xValue = std::round(allStrip[x].x);
			if (xValue >= 0 && xValue < WIDTH && yValue >= 0 && yValue < HEIGHT) {
				if (depthBuffer[yValue][xValue] < (1 / allStrip[x].depth)) {
					depthBuffer[yValue][xValue] = 1 / allStrip[x].depth;
					window.setPixelColour(xValue, yValue, col);
				}
			}
		}
	}
	//drawLine(window, triangle.v1(), triangle.v2(), colour);
}

void drawFilledTextureBottom(DrawingWindow& window, CanvasTriangle triangle, Colour colour, std::vector<std::vector<float>> &depthBuffer, TextureMap texture) {
	std::vector<glm::vec2> line1 = drawLine2(window, triangle.v0(), triangle.v1());
	std::vector<glm::vec2> line2 = drawLine2(window, triangle.v0(), triangle.v2());
	std::vector<glm::vec2> flat = drawLine2(window, triangle.v1(), triangle.v2());
	int triangleHeight = std::fabs(triangle.v2().y - triangle.v0().y);

	std::vector<float> line1correct;
	std::vector<float> line2correct;

	if (line1.size() == 0 || line2.size() == 0) return;

	/*line1correct = interpolateSingleFloats(line1[0][0], line1[line1.size() - 1][0], triangleHeight);
	line2correct = interpolateSingleFloats(line2[0][0], line2[line2.size() - 1][0], triangleHeight);

	std::vector<float> dLine1 = interpolateSingleFloats(triangle.v0().depth, triangle.v1().depth, triangleHeight);
	std::vector<float> dLine2 = interpolateSingleFloats(triangle.v0().depth, triangle.v2().depth, triangleHeight);*/

	std::vector<TexturePoint> tLine1 = interpolateTexture(triangle.v0().texturePoint, triangle.v1().texturePoint, triangleHeight + 1);
	std::vector<TexturePoint> tLine2 = interpolateTexture(triangle.v0().texturePoint, triangle.v2().texturePoint, triangleHeight + 1);

	std::vector<CanvasPoint> line1All = interpolateThreeElementValues(triangle.v0(), triangle.v1(), triangleHeight + 1);
	std::vector<CanvasPoint> line2All = interpolateThreeElementValues(triangle.v0(), triangle.v2(), triangleHeight + 1);

	//Texture perspective correction
	std::vector<CanvasPoint> farNear = furthestAndNearest(triangle);
	float z0 = farNear[0].depth;
	float z1 = farNear[1].depth;
	float c0 = farNear[0].texturePoint.y;
	float c1 = farNear[1].texturePoint.y;

	for (size_t y = 0; y < triangleHeight; y++) {
		int size = std::fabs(std::round((line1All[y].x - line2All[y].x))) * 2;
		/*std::vector<float> strip;
		strip = interpolateSingleFloats2(line1correct[y], line2correct[y], size);
		std::vector<float> dStrip;
		dStrip = interpolateSingleFloats(dLine1[y], dLine2[y], size);*/
		std::vector<TexturePoint> tStrip;
		tStrip = interpolateTexture(tLine1[y], tLine2[y], size);
		std::vector<CanvasPoint> allStrip = interpolateThreeElementValues(line1All[y], line2All[y], size);
		int q = y;

		for (size_t x = 0; x < allStrip.size(); x++) {
			int tIndex = texture.width * std::round(tStrip[x].y) + std::round(tStrip[x].x);
			//int c = (((c0 / z0) * (1 - q)) + ((c1 / z1) * q)) / (((1 / z0) * (1 - q)) + ((1 / z1) * q));
			//int tIndex = texture.width * c + std::round(tStrip[x].x);
			if (tIndex >= texture.pixels.size()) tIndex = texture.pixels.size() - 1;
			uint32_t col = texture.pixels[tIndex];

			int yValue = std::round(triangle.v0().y - y);
			int xValue = std::round(allStrip[x].x);
			if (xValue >= 0 && xValue < WIDTH && yValue >= 0 && yValue < HEIGHT) {
				if (depthBuffer[yValue][xValue] < (1 / allStrip[x].depth)) {
					depthBuffer[yValue][xValue] = 1 / allStrip[x].depth;
					window.setPixelColour(xValue, yValue, col);
				}
			}
		}
	}
}

void drawFilledTriangle(DrawingWindow& window, CanvasTriangle triangle, Colour colour, std::vector<std::vector<float>> &depthBuffer, TextureMap texture) {
	//window.clearPixels();
	bool yesTex = false;
	if (colour.name == "Cobbles" || colour.name == "Cobbles2") {
		yesTex = true;
		/*std::cout << "v0: " << triangle.v0().texturePoint << std::endl;
		std::cout << "v1: " << triangle.v1().texturePoint << std::endl;
		std::cout << "v2: " << triangle.v2().texturePoint << std::endl;*/
	}
	if (colour.name == "Glass") return;
	CanvasTriangle sorted = sort(triangle);
	//Colour white = Colour::Colour(255, 255, 255);
	float ratio = (sorted.v0().y - sorted.v2().y) / (sorted.v0().x - sorted.v2().x);
	float middleYDiff = std::fabs(sorted.v0().y - sorted.v1().y);
	float middleXDiff = std::fabs(sorted.v0().x - sorted.v1().x);
	CanvasPoint middle = CanvasPoint::CanvasPoint(std::round(sorted.v0().x + (middleYDiff * (1 / ratio))), sorted.v1().y);

	float dRatio = (middle.y - sorted.v0().y) / (sorted.v2().y - sorted.v0().y);
	middle.depth = sorted.v0().depth + dRatio * (sorted.v2().depth - sorted.v0().depth);

	//float tRatio = (sorted.v0().texturePoint.y - sorted.v2().texturePoint.y) / (sorted.v0().texturePoint.x - sorted.v2().texturePoint.x);
	//float textureX = std::round(sorted.v0().texturePoint.x + (middleYDiff * (1 / tRatio)));
	//float textureY = std::round(sorted.v0().texturePoint.y + (middleXDiff * (1 / tRatio)));

	float yRatio = (middle.y - sorted.v0().y) / (sorted.v2().y - sorted.v0().y);
	float textureX = sorted.v0().texturePoint.x + yRatio * (sorted.v2().texturePoint.x - sorted.v0().texturePoint.x);
	float xRatio = (middle.x - sorted.v0().x) / (sorted.v2().x - sorted.v0().x);
	float textureY = sorted.v0().texturePoint.y + xRatio * (sorted.v2().texturePoint.y - sorted.v0().texturePoint.y);
	middle.texturePoint = TexturePoint::TexturePoint(textureX, textureY);

	//Triangle 1
	CanvasTriangle triangle1 = CanvasTriangle::CanvasTriangle(sorted.v0(), sorted.v1(), middle);

	//Triangle 2
	CanvasTriangle triangle2 = CanvasTriangle::CanvasTriangle(sorted.v2(), middle, sorted.v1());

	//Right triangle with flat top
	if (sorted.v0().y == sorted.v1().y && (sorted.v0().x == sorted.v2().x || sorted.v1().x == sorted.v2().x)) {
		CanvasTriangle rTop = CanvasTriangle::CanvasTriangle(sorted.v2(), sorted.v0(), sorted.v1());
		if (!yesTex) drawFilledBottom(window, rTop, colour, depthBuffer);
		else drawFilledTextureBottom(window, rTop, colour, depthBuffer, texture);
	}

	//Right Triangle with flat bottom
	else if (sorted.v1().y == sorted.v2().y && (sorted.v1().x == sorted.v0().x || sorted.v2().x == sorted.v0().x)) {
		CanvasTriangle rBottom = CanvasTriangle::CanvasTriangle(sorted.v0(), sorted.v2(), sorted.v1());
		if (!yesTex) drawFilledTop(window, rBottom, colour, depthBuffer);
		else drawFilledTextureTop(window, rBottom, colour, depthBuffer, texture);
	}

	//Isoceles with flat top
	else if (sorted.v0().y == sorted.v1().y) {
		CanvasTriangle iTop = CanvasTriangle::CanvasTriangle(sorted.v2(), sorted.v0(), sorted.v1());
		if (!yesTex) drawFilledBottom(window, iTop, colour, depthBuffer);
		else drawFilledTextureBottom(window, iTop, colour, depthBuffer, texture);
	}

	//Isoceles with flat bottom
	else if (sorted.v1().y == sorted.v2().y) {
		if (!yesTex) drawFilledTop(window, sorted, colour, depthBuffer);
		else drawFilledTextureTop(window, sorted, colour, depthBuffer, texture);
	}

	//Other Triangles
	else {
		if (!yesTex) {
			drawFilledTop(window, triangle1, colour, depthBuffer);
			drawFilledBottom(window, triangle2, colour, depthBuffer);
		}
		else {
			drawFilledTextureTop(window, triangle1, colour, depthBuffer, texture);
			drawFilledTextureBottom(window, triangle2, colour, depthBuffer, texture);
		}
	}
	if (yesTex) {
		drawFilledTextureTop(window, triangle1, colour, depthBuffer, texture);
		drawFilledTextureBottom(window, triangle2, colour, depthBuffer, texture);
	}
	
	//drawTriangle(window, sorted, colour);
}

/*
void pointCloudRenderer(DrawingWindow& window, glm::vec3 cameraPosition, float focalLength, float imageScale, std::vector<ModelTriangle> models) {
	for (int i = 0; i < models.size(); i++) {
		float red = 255;
		float green = 255;
		float blue = 255;
		uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
		for (int j = 0; j < 3; j++) {
			glm::vec3 vectorPosition = models[i].vertices[j];
			CanvasPoint cPoint = getCanvasIntersectionPoint(cameraPosition, vectorPosition, focalLength, imageScale);
			window.setPixelColour(cPoint.x, cPoint.y, colour);
		}
	}
}
*/


void wireFrameRenderer(DrawingWindow& window, glm::vec3 cameraPosition, glm::mat3 cameraOrientation, float focalLength, float imageScale, std::vector<ModelTriangle> models) {
	window.clearPixels();
	for (int i = 0; i < models.size(); i++) {

		glm::vec3 vectorPosition1 = models[i].vertices[0];
		CanvasPoint cPoint1 = getCanvasIntersectionPoint(cameraPosition, cameraOrientation, vectorPosition1, focalLength, imageScale);
		glm::vec3 vectorPosition2 = models[i].vertices[1];
		CanvasPoint cPoint2 = getCanvasIntersectionPoint(cameraPosition, cameraOrientation, vectorPosition2, focalLength, imageScale);
		glm::vec3 vectorPosition3 = models[i].vertices[2];
		CanvasPoint cPoint3 = getCanvasIntersectionPoint(cameraPosition, cameraOrientation,vectorPosition3, focalLength, imageScale);
		CanvasTriangle triangle = CanvasTriangle::CanvasTriangle(cPoint1, cPoint2, cPoint3);
		drawTriangle(window, triangle, models[i].colour);
	}
}

void rasterisedRenderer(DrawingWindow& window, glm::vec3 cameraPosition, glm::mat3 cameraOrientation, float focalLength, float imageScale, std::vector<ModelTriangle> models, TextureMap texture) {
	window.clearPixels();
	std::vector<std::vector<float>> depthBuffer;
	for (int i = 0; i < HEIGHT; i++) {
		std::vector<float> temp;
		for (int j = 0; j < WIDTH; j++) {
			temp.push_back(0.0f);
		}
		depthBuffer.push_back(temp);
	}
	for (int i = 0; i < models.size(); i++) {
		//Colour white = Colour::Colour(255, 255, 255);

		std::vector<CanvasPoint> cPoints;
		for (int j = 0; j < 3; j++) {
			cPoints.push_back(getCanvasIntersectionPoint(cameraPosition, cameraOrientation, models[i].vertices[j], focalLength, imageScale));
		}
		CanvasTriangle triangle = CanvasTriangle::CanvasTriangle(cPoints[0], cPoints[1], cPoints[2]);
		triangle.v0().texturePoint = models[i].texturePoints[0];
		triangle.v1().texturePoint = models[i].texturePoints[1];
		triangle.v2().texturePoint = models[i].texturePoints[2];
		drawFilledTriangle(window, triangle, models[i].colour, depthBuffer, texture);
	}
}

glm::vec3 handleEvent(SDL_Event event, DrawingWindow &window, glm::vec3 cameraPosition, glm::mat3 cameraOrientation) {
	glm::vec3 newPosition = cameraPosition;
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) {
			std::cout << "LEFT" << std::endl;
		}
		else if (event.key.keysym.sym == SDLK_RIGHT) {
			std::cout << "RIGHT" << std::endl;
		}
		else if (event.key.keysym.sym == SDLK_UP) {
			std::cout << "UP" << std::endl;
		}
		else if (event.key.keysym.sym == SDLK_DOWN) {
			std::cout << "DOWN" << std::endl;
		}
		else if (event.key.keysym.sym == SDLK_f) {
			std::cout << "FORWARD" << std::endl;
		}
		else if (event.key.keysym.sym == SDLK_b) {
			std::cout << "BACKWARD" << std::endl;
		}
		newPosition = moveCamera(cameraPosition, cameraOrientation, event);
	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
	return newPosition;
}

glm::mat3 handleEvent2(SDL_Event event, DrawingWindow& window, glm::mat3 cameraOrientation) {
	glm::mat3 newOrientation = cameraOrientation;
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_w) {
			std::cout << "Rotate Up" << std::endl;
		}
		else if (event.key.keysym.sym == SDLK_s) {
			std::cout << "Rotate Down" << std::endl;
		}
		else if (event.key.keysym.sym == SDLK_d) {
			std::cout << "Rotate Right" << std::endl;
		}
		else if (event.key.keysym.sym == SDLK_a) {
			std::cout << "Rotate Left" << std::endl;
		}
		newOrientation = rotateCamera(cameraOrientation, event);
	}
	else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
	return newOrientation;
}

glm::vec3 handleEvent3(SDL_Event event, DrawingWindow& window, glm::vec3 lightSource) {
	glm::vec3 newPosition = lightSource;
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_j) {
			std::cout << "Light LEFT" << std::endl;
		}
		else if (event.key.keysym.sym == SDLK_l) {
			std::cout << "Light RIGHT" << std::endl;
		}
		else if (event.key.keysym.sym == SDLK_i) {
			std::cout << "Light UP" << std::endl;
		}
		else if (event.key.keysym.sym == SDLK_k) {
			std::cout << "Light DOWN" << std::endl;
		}
		else if (event.key.keysym.sym == SDLK_u) {
			std::cout << "Light FORWARD" << std::endl;
		}
		else if (event.key.keysym.sym == SDLK_o) {
			std::cout << "Light BACKWARD" << std::endl;
		}
		newPosition = moveLight(lightSource, event);
	}
	else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
	return newPosition;
}


std::vector<Colour> loadMaterials(const std::string& filename) {
	std::ifstream inputStream(filename);
	std::string nextLine;

	std::vector<Colour> colours;
	std::vector<glm::vec3> vectors;
	std::string cName;

	while (std::getline(inputStream, nextLine)) {
		if (nextLine.size() > 0 && nextLine.at(0) == 'n') {
			std::vector<std::string> line = split(nextLine, ' ');
			cName = line[1];
		}
		if (nextLine.size() > 0 && nextLine.at(0) == 'K') {
			std::vector<std::string> line = split(nextLine, ' ');
			float red = std::round(std::stof(line[1]) * 255);
			float green = std::round(std::stof(line[2]) * 255);
			float blue = std::round(std::stof(line[3]) * 255);
			Colour colour = Colour::Colour(cName, red, green, blue);
			colours.push_back(colour);
		}
	}

	inputStream.close();
	return colours;
}

Colour findColour(std::vector<Colour> pallete, std::string name) {
	Colour colour = Colour::Colour(255, 255, 255);
	for (int i = 0; i < pallete.size(); i++) {
		if (name == pallete[i].name) return pallete[i];
	}
	return colour;
}

std::vector<ModelTriangle> loadModel(const std::string& filename, float scallingFactor, std::vector<Colour> palette, TextureMap textureMap) {
	std::ifstream inputStream(filename);
	std::string nextLine;

	std::vector<ModelTriangle> models;
	std::vector<glm::vec3> vectors;
	std::vector<TexturePoint> textures;
	Colour colour;

	while (std::getline(inputStream, nextLine)) {
		if (nextLine.size() > 0 && nextLine.at(0) == 'u') {
			std::vector<std::string> line = split(nextLine, ' ');
			std::string cName = line[1];
			colour = findColour(palette, cName);
		}
		if (nextLine.size() > 0 && nextLine.at(0) == 'v' && nextLine.at(1) != 't' ){
			std::vector<std::string> line = split(nextLine, ' ');
			float point1 = std::stof(line[1]) * scallingFactor;
			float point2 = std::stof(line[2]) * scallingFactor;
			float point3 = std::stof(line[3]) * scallingFactor;
			glm::vec3 vector = glm::vec3(point1, point2, point3);
			vectors.push_back(vector);
		}
		if (nextLine.size() > 0 && nextLine.at(0) == 'v' && nextLine.at(1) == 't') {
			std::vector<std::string> line = split(nextLine, ' ');
			float point1 = std::stof(line[1]);
			float point2 = std::stof(line[2]);
			TexturePoint texture = TexturePoint::TexturePoint(std::round(point1 * textureMap.width), std::round(point2 * textureMap.height));
			textures.push_back(texture);
		}
		if (nextLine.size() > 0 && nextLine.at(0) == 'f') {
			std::vector<std::string> line = split(nextLine, ' ');
			std::vector<std::string> part1 = split(line[1], '/');
			std::vector<std::string> part2 = split(line[2], '/');
			std::vector<std::string> part3 = split(line[3], '/');

			float point1 = std::stof(part1[0]);
			float point2 = std::stof(part2[0]);
			float point3 = std::stof(part3[0]);

			float tex1;
			float tex2;
			float tex3;
			TexturePoint t0;
			TexturePoint t1;
			TexturePoint t2;
			bool yesTex = false;
			//if (part1.size() == 2 && part2.size() == 2 && part3.size() == 2) yesTex = true;
			//if (part1[1].at(0) != ' ' && part2[1].at(0) != ' ' && part3[1].at(0) != ' ') yesTex = true;
			if (colour.name == "Cobbles" || colour.name == "Cobbles2") yesTex = true;
			if (yesTex) {
				tex1 = std::stof(part1[1]);
				tex2 = std::stof(part2[1]);
				tex3 = std::stof(part3[1]);
				t0 = textures[tex1 - 1];
				t1 = textures[tex2 - 1];
				t2 = textures[tex3 - 1];
			}
			glm::vec3 v0 = vectors[point1 - 1];
			glm::vec3 v1 = vectors[point2 - 1];
			glm::vec3 v2 = vectors[point3 - 1];
			//std::cout << colour << std::endl;
			ModelTriangle model = ModelTriangle::ModelTriangle(v0, v1, v2, colour);
			
			model.normal = glm::cross((v1 - v0), (v2 - v0));
			//model.normal = glm::cross((v2 - v0), (v1 - v0));
			
			if (yesTex) {
				model.texturePoints = std::array<TexturePoint, 3>{t0, t1, t2};
			}
			models.push_back(model);
		}

	}

	inputStream.close();
	return models;
}

float min_el(std::vector<float> v) {
	float min = v[0];
	for (int i = 1; i < v.size(); i++) {
		if (v[i] < min) min = v[i];
	}
	return min;
}

float max_el(std::vector<float> v) {
	float max = v[0];
	for (int i = 1; i < v.size(); i++) {
		if (v[i] > max) max = v[i];
	}
	return max;
}

std::vector<glm::vec3> modelDimensions(std::vector<ModelTriangle> models) {
	std::vector<float> xs;
	std::vector<float> ys;
	std::vector<float> zs;
	for (int i = 0; i < models.size(); i++) {
		for (int j = 0; j < 3; j++) {
			xs.push_back(models[i].vertices[j].x);
			ys.push_back(models[i].vertices[j].y);
			zs.push_back(models[i].vertices[j].z);
		}
	}
	float xmin = min_el(xs);
	//std::cout << "min x: " << xmin << std::endl;
	float xmax = max_el(xs);
	//std::cout << "max x: " << xmax << std::endl;
	float ymin = min_el(ys);
	//std::cout << "min y: " << ymin << std::endl;
	float ymax = max_el(ys);
	//std::cout << "max y: " << ymax << std::endl;
	float zmin = min_el(zs);
	//std::cout << "min z: " << zmin << std::endl;
	float zmax = max_el(zs);
	//std::cout << "max z: " << zmax << std::endl;
	float width = std::fabs(xmax - xmin);
	float height = std::fabs(ymax - ymin);
	float depth = std::fabs(zmax - zmin);
	std::vector<glm::vec3> dimensions;
	glm::vec3 sizes = glm::vec3(width, height, depth);
	glm::vec3 mins = glm::vec3(xmin, ymin, zmin);
	dimensions.push_back(sizes);
	dimensions.push_back(mins);
	return dimensions;
}

//photon firing function
glm::vec3 photonTrace(glm::vec3 start, glm::vec3 rayDirection, float intensity, std::vector<ModelTriangle> models) {
	RayTriangleIntersection  intersection = getClosestIntersection(start, rayDirection, models);
	glm::vec3 reflectionVector = glm::normalize(rayDirection) - 2.0f * glm::normalize(intersection.intersectedTriangle.normal) * glm::dot(glm::normalize(rayDirection), glm::normalize(intersection.intersectedTriangle.normal));
	if (intersection.intersectedTriangle.colour.name == "Mirror") {
		intersection = getClosestReflection(intersection.intersectionPoint, reflectionVector, models);
	}
	if (intersection.intersectedTriangle.colour.name == "Glass") {
		glm::vec3 normal = intersection.intersectedTriangle.normal;
		float glassRI = 1.5;
		float airRI = 1.0;
		float N1 = airRI;
		float N2 = glassRI;
		float incidentAngle = glm::dot(normal,rayDirection);
		float RI = N1 / N2;
		float theta1 = incidentAngle;
		glm::vec3 I = rayDirection;
		float c1 = glm::dot(normal, I);
		float refractionAngle = std::asin(std::sin(theta1) * RI);
		float theta2 = refractionAngle;
		float c2 = (float)std::sqrt(1 - (std::pow(RI, 2) * (1 - std::pow(std::cos(theta1), 2))));
		glm::vec3 tRay = RI * I + (RI * c1 - c2) * normal;
		glm::vec3 bias = normal * 0.0f;
		intersection = getClosestRefraction(intersection.intersectionPoint - bias, tRay, models);
	}
	float fabsorbProb = containAngle(glm::dot(intersection.intersectedTriangle.normal, rayDirection));
	float p = (rand() % 100 + 1) / 100.0;
	if (p > fabsorbProb) {
		return photonTrace(intersection.intersectionPoint, glm::normalize(reflectionVector), intensity, models);
	}
	else return intersection.intersectionPoint;
}

std::vector<std::vector<std::vector<float>>> photonMap(glm::vec3 lightSource, std::vector<glm::vec3> dimensions, float photons, float intensity, int aggRadius, int scale, std::vector<ModelTriangle> models) {
	float width = dimensions[0][0] * scale +1;
	float height = dimensions[0][1] * scale +1;
	float depth = dimensions[0][2] * scale +1;
	float xmin = dimensions[1][0];
	float ymin = dimensions[1][1];
	float zmin = dimensions[1][2];
	std::vector<std::vector<std::vector<float>>> pMap;
	for (int i = 0; i <= width; i++) {
		std::vector<std::vector<float>> two;
		for (int j = 0; j <= height; j++) {
			std::vector<float> one;
			for (int k = 0; k <= depth; k++) {
				one.push_back(0);
			}
			two.push_back(one);
		}
		pMap.push_back(two);
	}
	//ray for loop
	for (int i = 0; i < photons; i++) {
		float randX = rand() % 100;
		float randY = rand() % 100;
		float randZ = rand() % 100;
		glm::vec3 lightRay = glm::normalize(lightSource - glm::vec3(randX, randY, randZ));
		glm::vec3 fabsorbPoint = photonTrace(lightSource, lightRay, intensity, models);
		int xPos = std::round((fabsorbPoint[0] -xmin) * scale);
		int yPos = std::round((fabsorbPoint[1] -ymin) * scale);
		int zPos = std::round((fabsorbPoint[2] -zmin) * scale);
		pMap[xPos][yPos][zPos] += intensity;
	}

	//aggregation 
	std::vector<std::vector<std::vector<float>>> agg = pMap;
	for (int i = 0; i <= width; i++) {
		for (int j = 0; j <= height; j++) {
			for (int k = 0; k <= depth; k++) {
				float average = 0;
				float count = 0;
				for (int b = -aggRadius; b <= aggRadius; b++) {
					for (int n = -aggRadius; n <= aggRadius; n++) {
						for (int m = -aggRadius; m <= aggRadius; m++) {
							int x = i + b;
							int y = j + n;
							int z = k + m;
							if (x >= 0 && x <= width && y >= 0 && y <= height && z >= 0 && z <= depth && std::fabs(b) + std::fabs(n) + std::fabs(m) <= aggRadius) {
								count += 1;
								average += pMap[x][y][z];
							}
						}
					}
				}
				agg[i][j][k] = average / count;
			}
		}
	}
	pMap = agg;
	return pMap;
}

//lightSource = moveLight(lightSource, event);
//cameraPosition = moveCamera(cameraPosition, cameraOrientation, event);
//cameraOrientation = rotateCamera(cameraOrientation, event);
//SDL_Event event;
//event.type == SDL_KEYDOWN
//rasterisedRenderer(window, cameraPosition, cameraOrientation, focalLength, imageScale, models, texture);
//rayTrace(window, cameraPosition, cameraOrientation, lightSource, focalLength, imageScale, shadingCode, shadowCode, models);
//wireFrameRenderer(window, cameraPosition, cameraOrientation, focalLength, imageScale, models);
//window.renderFrame();

void autoViewing1(std::vector<ModelTriangle> models, TextureMap texture, DrawingWindow& window, float focalLength, float imageScale, SDL_Event event) {
	glm::vec3 cameraPosition = glm::vec3(0.0, 0.0, 4.0);
	glm::mat3 cameraOrientation = glm::mat3(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
	glm::vec3 lightSource = glm::vec3(-0.1, .3, 0.1);
	int count = 0;

	wireFrameRenderer(window, cameraPosition, cameraOrientation, focalLength, imageScale, models);
	for (int i = 0; i < 15; i++) {
		std::this_thread::sleep_for(std::chrono::milliseconds(50));
		window.renderFrame();
		count += 1;
		std::string output = "../../../output/scene1/pic" + std::to_string(count) + std::string(".bmp");
		window.saveBMP(output);
	}
	
	//Rotate left 
	event.key.keysym.sym = SDLK_d;
	for (int i = 0; i < 15; i++) {
		cameraPosition = moveCamera(cameraPosition, cameraOrientation, event);
		cameraOrientation = rotateCamera(cameraOrientation, event);
		rasterisedRenderer(window, cameraPosition, cameraOrientation, focalLength, imageScale, models, texture);
		window.renderFrame();
		count += 1;
		std::string output = "../../../output/scene1/pic" + std::to_string(count) + std::string(".bmp");
		window.saveBMP(output);
	}
	//Rotate Up
	event.key.keysym.sym = SDLK_s;
	for (int i = 0; i < 15; i++) {
		cameraPosition = moveCamera(cameraPosition, cameraOrientation, event);
		cameraOrientation = rotateCamera(cameraOrientation, event);
		rasterisedRenderer(window, cameraPosition, cameraOrientation, focalLength, imageScale, models, texture);
		window.renderFrame();
		count += 1;
		std::string output = "../../../output/scene1/pic" + std::to_string(count) + std::string(".bmp");
		window.saveBMP(output);
	}

	rayTrace(window, cameraPosition, cameraOrientation, lightSource, focalLength, imageScale, 1, 0, models);
	std::cout << "Flat Shading" << std::endl;
	for (int i = 0; i < 60; i++) {
		//std::this_thread::sleep_for(std::chrono::seconds(1));
		window.renderFrame();
		count += 1;
		std::string output = "../../../output/scene1/pic" + std::to_string(count) + std::string(".bmp");
		window.saveBMP(output);
	}

	rayTrace(window, cameraPosition, cameraOrientation, lightSource, focalLength, imageScale, 2, 0, models);
	std::cout << "Goraud Shading" << std::endl;
	for (int i = 0; i < 75; i++) {
		//std::this_thread::sleep_for(std::chrono::seconds(1));
		window.renderFrame();
		count += 1;
		std::string output = "../../../output/scene1/pic" + std::to_string(count) + std::string(".bmp");
		window.saveBMP(output);
	}

	rayTrace(window, cameraPosition, cameraOrientation, lightSource, focalLength, imageScale, 3, 0, models);
	std::cout << "Phong" << std::endl;
	for (int i = 0; i < 75; i++) {
		std::this_thread::sleep_for(std::chrono::milliseconds(500));
		window.renderFrame();
		count += 1;
		std::string output = "../../../output/scene1/pic" + std::to_string(count) + std::string(".bmp");
		window.saveBMP(output);
	}

	//total frames is 255
}

void autoViewing2(std::vector<ModelTriangle> models, TextureMap texture, DrawingWindow& window, float focalLength, float imageScale, SDL_Event event) {
	glm::vec3 cameraPosition = glm::vec3(0.0, 0.0, 3.75);
	glm::mat3 cameraOrientation = glm::mat3(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
	glm::vec3 lightSource = glm::vec3(0.0, .3, 0.0);
	int count = 0;

	wireFrameRenderer(window, cameraPosition, cameraOrientation, focalLength, imageScale, models);
	for (int i = 0; i < 15; i++) {
		std::this_thread::sleep_for(std::chrono::milliseconds(50));
		window.renderFrame();
		count += 1;
		std::string output = "../../../output/scene2/pic" + std::to_string(count) + std::string(".bmp");
		window.saveBMP(output);
	}
	
	rasterisedRenderer(window, cameraPosition, cameraOrientation, focalLength, imageScale, models, texture);
	for (int i = 0; i < 30; i++) {
		//std::this_thread::sleep_for(std::chrono::milliseconds(500));
		window.renderFrame();
		count += 1;
		std::string output = "../../../output/scene2/pic" + std::to_string(count) + std::string(".bmp");
		window.saveBMP(output);
	}

	rayTrace(window, cameraPosition, cameraOrientation, lightSource, focalLength, imageScale, 1, 1, models);
	for (int i = 0; i < 60; i++) {
		//std::this_thread::sleep_for(std::chrono::milliseconds(500));
		window.renderFrame();
		count += 1;
		std::string output = "../../../output/scene2/pic" + std::to_string(count) + std::string(".bmp");
		window.saveBMP(output);
	}

	//Rotate left 
	event.key.keysym.sym = SDLK_d;
	for (int i = 0; i < 15; i++) {
		cameraPosition = moveCamera(cameraPosition, cameraOrientation, event);
		cameraOrientation = rotateCamera(cameraOrientation, event);
		rasterisedRenderer(window, cameraPosition, cameraOrientation, focalLength, imageScale, models, texture);
		window.renderFrame();
		count += 1;
		std::string output = "../../../output/scene2/pic" + std::to_string(count) + std::string(".bmp");
		window.saveBMP(output);
	}

	rayTrace(window, cameraPosition, cameraOrientation, lightSource, focalLength, imageScale, 1, 1, models);
	for (int i = 0; i < 60; i++) {
		//std::this_thread::sleep_for(std::chrono::milliseconds(500));
		window.renderFrame();
		count += 1;
		std::string output = "../../../output/scene2/pic" + std::to_string(count) + std::string(".bmp");
		window.saveBMP(output);
	}
	//Total Frames 180
}

void autoViewing3(std::vector<ModelTriangle> models, TextureMap texture, DrawingWindow& window, float focalLength, float imageScale, SDL_Event event) {
	glm::vec3 cameraPosition = glm::vec3(0.0, 0.0, 4.0);
	glm::mat3 cameraOrientation = glm::mat3(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
	glm::vec3 lightSource = glm::vec3(0.0, .3, 0.0);
	int count = 0;

	wireFrameRenderer(window, cameraPosition, cameraOrientation, focalLength, imageScale, models);
	for (int i = 0; i < 15; i++) {
		std::this_thread::sleep_for(std::chrono::milliseconds(50));
		window.renderFrame();
		count += 1;
		std::string output = "../../../output/scene3/pic" + std::to_string(count) + std::string(".bmp");
		window.saveBMP(output);
	}

	//Rotate Up
	event.key.keysym.sym = SDLK_s;
	for (int i = 0; i < 15; i++) {
		cameraPosition = moveCamera(cameraPosition, cameraOrientation, event);
		cameraOrientation = rotateCamera(cameraOrientation, event);
		rasterisedRenderer(window, cameraPosition, cameraOrientation, focalLength, imageScale, models, texture);
		window.renderFrame();
		count += 1;
		std::string output = "../../../output/scene3/pic" + std::to_string(count) + std::string(".bmp");
		window.saveBMP(output);
	}

	rayTrace(window, cameraPosition, cameraOrientation, lightSource, focalLength, imageScale, 1, 0, models);
	std::cout << "No Shadows" << std::endl;
	for (int i = 0; i < 75; i++) {
		//std::this_thread::sleep_for(std::chrono::seconds(1));
		window.renderFrame();
		count += 1;
		std::string output = "../../../output/scene3/pic" + std::to_string(count) + std::string(".bmp");
		window.saveBMP(output);
	}

	rayTrace(window, cameraPosition, cameraOrientation, lightSource, focalLength, imageScale, 1, 1, models);
	std::cout << "Hard Shadows" << std::endl;
	for (int i = 0; i < 75; i++) {
		//std::this_thread::sleep_for(std::chrono::seconds(1));
		window.renderFrame();
		count += 1;
		std::string output = "../../../output/scene3/pic" + std::to_string(count) + std::string(".bmp");
		window.saveBMP(output);
	}

	rayTrace(window, cameraPosition, cameraOrientation, lightSource, focalLength, imageScale, 1, 2, models);
	std::cout << "Soft Shadows" << std::endl;
	for (int i = 0; i < 75; i++) {
		std::this_thread::sleep_for(std::chrono::milliseconds(500));
		window.renderFrame();
		count += 1;
		std::string output = "../../../output/scene3/pic" + std::to_string(count) + std::string(".bmp");
		window.saveBMP(output);
	}
	//245 Total Frames
}

void autoViewing4(std::vector<ModelTriangle> models, TextureMap texture, DrawingWindow& window, float focalLength, float imageScale, SDL_Event event) {
	glm::vec3 cameraPosition = glm::vec3(0.0, 0.0, 4.0);
	glm::mat3 cameraOrientation = glm::mat3(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
	glm::vec3 lightSource = glm::vec3(0.0, .05, 0.0);
	int count = 0;

	wireFrameRenderer(window, cameraPosition, cameraOrientation, focalLength, imageScale, models);
	for (int i = 0; i < 15; i++) {
		std::this_thread::sleep_for(std::chrono::milliseconds(50));
		window.renderFrame();
	}

	//Rotate Up
	event.key.keysym.sym = SDLK_s;
	for (int i = 0; i < 15; i++) {
		cameraPosition = moveCamera(cameraPosition, cameraOrientation, event);
		cameraOrientation = rotateCamera(cameraOrientation, event);
		rasterisedRenderer(window, cameraPosition, cameraOrientation, focalLength, imageScale, models, texture);
		window.renderFrame();
		count += 1;
		std::string output = "../../../output/scene4/pic" + std::to_string(count) + std::string(".bmp");
		window.saveBMP(output);
	}

	//Rotate left 
	event.key.keysym.sym = SDLK_a;
	for (int i = 0; i < 30; i++) {
		cameraPosition = moveCamera(cameraPosition, cameraOrientation, event);
		cameraOrientation = rotateCamera(cameraOrientation, event);
		rasterisedRenderer(window, cameraPosition, cameraOrientation, focalLength, imageScale, models, texture);
		window.renderFrame();
		count += 1;
		std::string output = "../../../output/scene4/pic" + std::to_string(count) + std::string(".bmp");
		window.saveBMP(output);
	}
	cameraPosition = glm::vec3(0.0, 0.0, 4.0);
	cameraOrientation = glm::mat3(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);

	//Rotate Up
	event.key.keysym.sym = SDLK_s;
	for (int i = 0; i < 30; i++) {
		cameraPosition = moveCamera(cameraPosition, cameraOrientation, event);
		cameraOrientation = rotateCamera(cameraOrientation, event);
		rasterisedRenderer(window, cameraPosition, cameraOrientation, focalLength, imageScale, models, texture);
		window.renderFrame();
		count += 1;
		std::string output = "../../../output/scene4/pic" + std::to_string(count) + std::string(".bmp");
		window.saveBMP(output);
	}

	rayTrace(window, cameraPosition, cameraOrientation, lightSource, focalLength, imageScale, 1, 0, models);
	std::cout << "Flat Shading" << std::endl;
	for (int i = 0; i < 75; i++) {
		//std::this_thread::sleep_for(std::chrono::seconds(1));
		window.renderFrame();
		count += 1;
		std::string output = "../../../output/scene4/pic" + std::to_string(count) + std::string(".bmp");
		window.saveBMP(output);
	}

	rayTrace(window, cameraPosition, cameraOrientation, lightSource, focalLength, imageScale, 2, 0, models);
	std::cout << "Goraud Shading" << std::endl;
	for (int i = 0; i < 75; i++) {
		//std::this_thread::sleep_for(std::chrono::seconds(1));
		window.renderFrame();
		count += 1;
		std::string output = "../../../output/scene4/pic" + std::to_string(count) + std::string(".bmp");
		window.saveBMP(output);
	}

	rayTrace(window, cameraPosition, cameraOrientation, lightSource, focalLength, imageScale, 3, 0, models);
	std::cout << "Phong" << std::endl;
	for (int i = 0; i < 75; i++) {
		std::this_thread::sleep_for(std::chrono::milliseconds(500));
		window.renderFrame();
		count += 1;
		std::string output = "../../../output/scene4/pic" + std::to_string(count) + std::string(".bmp");
		window.saveBMP(output);
	}
	// 300 Total Frames 
} 

void autoViewing5(std::vector<ModelTriangle> models, TextureMap texture, DrawingWindow& window, float focalLength, float imageScale, SDL_Event event) {
	glm::vec3 cameraPosition = glm::vec3(0.0, 0.0, 4.0);
	glm::mat3 cameraOrientation = glm::mat3(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
	glm::vec3 lightSource = glm::vec3(0.0, .3, 0.0);
	int count = 0;
	
	wireFrameRenderer(window, cameraPosition, cameraOrientation, focalLength, imageScale, models);
	for (int i = 0; i < 15; i++) {
		//std::this_thread::sleep_for(std::chrono::milliseconds(50));
		window.renderFrame();
		count += 1;
		std::string output = "../../../output/scene5/pic" + std::to_string(count) + std::string(".bmp");
		window.saveBMP(output);
	}
	
	rasterisedRenderer(window, cameraPosition, cameraOrientation, focalLength, imageScale, models, texture);
	for (int i = 0; i < 30; i++) {
		//std::this_thread::sleep_for(std::chrono::milliseconds(500));
		window.renderFrame();
		count += 1;
		std::string output = "../../../output/scene5/pic" + std::to_string(count) + std::string(".bmp");
		window.saveBMP(output);
	}

	//without Photon Map
	rayTrace(window, cameraPosition, cameraOrientation, lightSource, focalLength, imageScale, 3, 2, models);
	for (int i = 0; i < 60; i++) {
		//std::this_thread::sleep_for(std::chrono::seconds(1));
		window.renderFrame();
		count += 1;
		std::string output = "../../../output/scene5/pic" + std::to_string(count) + std::string(".bmp");
		window.saveBMP(output);
	}

	//with Photon Map
	int scale = 100;
	std::vector<glm::vec3> dimensions = modelDimensions(models);
	std::vector<std::vector<std::vector<float>>> pMap = photonMap(lightSource, dimensions, 10000, 1, 7, scale, models);
	std::cout << "Photon Mapping" << std::endl;
	rayTraceP(window, cameraPosition, cameraOrientation, lightSource, focalLength, imageScale, 3, 2, models, scale, dimensions, pMap);
	for (int i = 0; i < 60; i++) {
		std::this_thread::sleep_for(std::chrono::seconds(1));
		window.renderFrame();
		count += 1;
		std::string output = "../../../output/scene5/pic" + std::to_string(count) + std::string(".bmp");
		window.saveBMP(output);
	}

	//240 Total Frames
}

void autoViewing6(std::vector<ModelTriangle> models, TextureMap texture, DrawingWindow& window, float focalLength, float imageScale, SDL_Event event) {
	glm::vec3 cameraPosition = glm::vec3(0.0, 0.0, 4.0);
	glm::mat3 cameraOrientation = glm::mat3(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
	glm::vec3 lightSource = glm::vec3(0.0, .3, 0.1);
	int count = 0;

	wireFrameRenderer(window, cameraPosition, cameraOrientation, focalLength, imageScale, models);
	for (int i = 0; i < 30; i++) {
		//std::this_thread::sleep_for(std::chrono::milliseconds(50));
		window.renderFrame();
		count += 1;
		std::string output = "../../../output/scene6/pic" + std::to_string(count) + std::string(".bmp");
		window.saveBMP(output);
	}

	//Move Right
	event.key.keysym.sym = SDLK_RIGHT;
	for (int i = 0; i < 30; i++) {
		cameraPosition = moveCamera(cameraPosition, cameraOrientation, event);
		rasterisedRenderer(window, cameraPosition, cameraOrientation, focalLength, imageScale, models, texture);
		window.renderFrame();
		count += 1;
		std::string output = "../../../output/scene6/pic" + std::to_string(count) + std::string(".bmp");
		window.saveBMP(output);
	}

	//Move Left
	event.key.keysym.sym = SDLK_LEFT;
	for (int i = 0; i < 30; i++) {
		cameraPosition = moveCamera(cameraPosition, cameraOrientation, event);
		rasterisedRenderer(window, cameraPosition, cameraOrientation, focalLength, imageScale, models, texture);
		window.renderFrame();
		count += 1;
		std::string output = "../../../output/scene6/pic" + std::to_string(count) + std::string(".bmp");
		window.saveBMP(output);
	}

	//Move UP
	event.key.keysym.sym = SDLK_UP;
	for (int i = 0; i < 30; i++) {
		cameraPosition = moveCamera(cameraPosition, cameraOrientation, event);
		rasterisedRenderer(window, cameraPosition, cameraOrientation, focalLength, imageScale, models, texture);
		window.renderFrame();
		count += 1;
		std::string output = "../../../output/scene6/pic" + std::to_string(count) + std::string(".bmp");
		window.saveBMP(output);
	}

	//Move Down
	event.key.keysym.sym = SDLK_DOWN;
	for (int i = 0; i < 30; i++) {
		cameraPosition = moveCamera(cameraPosition, cameraOrientation, event);
		rasterisedRenderer(window, cameraPosition, cameraOrientation, focalLength, imageScale, models, texture);
		window.renderFrame();
		count += 1;
		std::string output = "../../../output/scene6/pic" + std::to_string(count) + std::string(".bmp");
		window.saveBMP(output);
	}

	//Move Back
	event.key.keysym.sym = SDLK_b;
	for (int i = 0; i < 30; i++) {
		cameraPosition = moveCamera(cameraPosition, cameraOrientation, event);
		rasterisedRenderer(window, cameraPosition, cameraOrientation, focalLength, imageScale, models, texture);
		window.renderFrame();
		count += 1;
		std::string output = "../../../output/scene6/pic" + std::to_string(count) + std::string(".bmp");
		window.saveBMP(output);
	}

	//Move Forward
	event.key.keysym.sym = SDLK_f;
	for (int i = 0; i < 30; i++) {
		cameraPosition = moveCamera(cameraPosition, cameraOrientation, event);
		rasterisedRenderer(window, cameraPosition, cameraOrientation, focalLength, imageScale, models, texture);
		window.renderFrame();
		count += 1;
		std::string output = "../../../output/scene6/pic" + std::to_string(count) + std::string(".bmp");
		window.saveBMP(output);
	}

	//Rotate left 
	event.key.keysym.sym = SDLK_d;
	for (int i = 0; i < 15; i++) {
		cameraPosition = moveCamera(cameraPosition, cameraOrientation, event);
		cameraOrientation = rotateCamera(cameraOrientation, event);
		rasterisedRenderer(window, cameraPosition, cameraOrientation, focalLength, imageScale, models, texture);
		window.renderFrame();
		count += 1;
		std::string output = "../../../output/scene6/pic" + std::to_string(count) + std::string(".bmp");
		window.saveBMP(output);
	}

	//Rotate right
	event.key.keysym.sym = SDLK_a;
	for (int i = 0; i < 15; i++) {
		cameraPosition = moveCamera(cameraPosition, cameraOrientation, event);
		cameraOrientation = rotateCamera(cameraOrientation, event);
		rasterisedRenderer(window, cameraPosition, cameraOrientation, focalLength, imageScale, models, texture);
		window.renderFrame();
		count += 1;
		std::string output = "../../../output/scene6/pic" + std::to_string(count) + std::string(".bmp");
		window.saveBMP(output);
	}

	//Rotate Down
	event.key.keysym.sym = SDLK_w;
	for (int i = 0; i < 15; i++) {
		cameraPosition = moveCamera(cameraPosition, cameraOrientation, event);
		cameraOrientation = rotateCamera(cameraOrientation, event);
		rasterisedRenderer(window, cameraPosition, cameraOrientation, focalLength, imageScale, models, texture);
		window.renderFrame();
		count += 1;
		std::string output = "../../../output/scene6/pic" + std::to_string(count) + std::string(".bmp");
		window.saveBMP(output);
	}

	//Rotate Up
	event.key.keysym.sym = SDLK_s;
	for (int i = 0; i < 15; i++) {
		cameraPosition = moveCamera(cameraPosition, cameraOrientation, event);
		cameraOrientation = rotateCamera(cameraOrientation, event);
		rasterisedRenderer(window, cameraPosition, cameraOrientation, focalLength, imageScale, models, texture);
		window.renderFrame();
		count += 1;
		std::string output = "../../../output/scene6/pic" + std::to_string(count) + std::string(".bmp");
		window.saveBMP(output);
	}

	rasterisedRenderer(window, cameraPosition, cameraOrientation, focalLength, imageScale, models, texture);
	for (int i = 0; i < 30; i++) {
		//std::this_thread::sleep_for(std::chrono::milliseconds(500));
		window.renderFrame();
		count += 1;
		std::string output = "../../../output/scene6/pic" + std::to_string(count) + std::string(".bmp");
		window.saveBMP(output);
	}

	rayTrace(window, cameraPosition, cameraOrientation, lightSource, focalLength, imageScale, 1, 1, models);
	std::cout << "Flat Shading" << std::endl;
	for (int i = 0; i < 60; i++) {
		//std::this_thread::sleep_for(std::chrono::seconds(1));
		window.renderFrame();
		count += 1;
		std::string output = "../../../output/scene6/pic" + std::to_string(count) + std::string(".bmp");
		window.saveBMP(output);
	}

	rayTrace(window, cameraPosition, cameraOrientation, lightSource, focalLength, imageScale, 2, 1, models);
	std::cout << "Goraud Shading" << std::endl;
	for (int i = 0; i < 75; i++) {
		//std::this_thread::sleep_for(std::chrono::seconds(1));
		window.renderFrame();
		count += 1;
		std::string output = "../../../output/scene6/pic" + std::to_string(count) + std::string(".bmp");
		window.saveBMP(output);
	}

	rayTrace(window, cameraPosition, cameraOrientation, lightSource, focalLength, imageScale, 3, 1, models);
	std::cout << "Phong" << std::endl;
	for (int i = 0; i < 75; i++) {
		std::this_thread::sleep_for(std::chrono::milliseconds(500));
		window.renderFrame();
		count += 1;
		std::string output = "../../../output/scene6/pic" + std::to_string(count) + std::string(".bmp");
		window.saveBMP(output);
	}

	//total frames is 255
}


void autoViewing(std::vector<std::vector<ModelTriangle>> models, std::vector<TextureMap> textures, DrawingWindow& window, float focalLength, float imageScale, SDL_Event event) {
	//autoViewing1(models[0], textures[0], window, focalLength, imageScale, event);
	autoViewing6(models[0], textures[0], window, focalLength, imageScale, event); //Actually scene 1 as original scene 1 isn't used. 
	//autoViewing2(models[1], textures[0], window, focalLength, imageScale, event);
	//autoViewing3(models[2], textures[0], window, focalLength, imageScale, event);
	//autoViewing4(models[3], textures[0], window, focalLength, imageScale, event);
	//autoViewing5(models[4], textures[0], window, focalLength, imageScale, event);
}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;

	std::vector<TextureMap> textures;
	textures.push_back(TextureMap::TextureMap("../../../src/texture.ppm"));

	std::vector<Colour> materials;
	materials = loadMaterials("../../../src/scene.mtl");

	std::vector<std::vector<ModelTriangle>> models;
	Colour white = Colour::Colour(255, 255, 255);
	models.push_back(loadModel("../../../src/textured-cornell-box.obj", .25, materials, textures[0]));
	models.push_back(loadModel("../../../src/IceScene.obj", .20, materials, textures[0]));
	models.push_back(loadModel("../../../src/EarthScene.obj", .20, materials, textures[0]));
	models.push_back(loadModel("../../../src/AirScene.obj", .25, materials, textures[0]));
	models.push_back(loadModel("../../../src/CBUltimate2.obj", .25, materials, textures[0]));
	

	glm::vec3 cameraPosition = glm::vec3(0.0, 0.0, 4.0);
	float focalLength = 2.0;
	float imageScale = 420;

	glm::mat3 cameraOrientation = glm::mat3(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);

	glm::vec3 rayDirection = glm::vec3(-0.1, -0.1, -2.0);
	//RayTriangleIntersection test = getClosestIntersection(cameraPosition, rayDirection, models);
	
	glm::vec3 lightSource = glm::vec3(0.0, .3, 0.0);
	
	int renderingCode = 2; // 0 = no rendering, 1 = wire frame, 2 = rasterized, 3 = ray trace
	int shadingCode = 1; // 1 = Flat, 2 = Gouraud, 3 = phong
	int shadowCode = 1; // 0 = No Shadows, 1 = Hard, 2 = Soft
	int sceneCode = 5; //Switch between different loaded scenes
	bool render = true;

	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) {
			if (event.type == SDL_KEYDOWN) render = true; // only render if there has been an update
			//Rendering Switch
			if (event.type == SDL_KEYDOWN && event.key.keysym.sym == SDLK_y && renderingCode != 1) {
				std::cout << "switching to wireframe" << std::endl;
				renderingCode = 1;
			}
			else if (event.type == SDL_KEYDOWN && event.key.keysym.sym == SDLK_r && renderingCode != 2) {
				std::cout << "switching to rasterized" << std::endl;
				renderingCode = 2;
			}
			else if (event.type == SDL_KEYDOWN && event.key.keysym.sym == SDLK_t && renderingCode != 3) {
				std::cout << "switching to ray traced" << std::endl;
				renderingCode = 3;
			}
			//Shading Switch
			if (event.type == SDL_KEYDOWN && event.key.keysym.sym == SDLK_c && shadingCode != 1) {
				std::cout << "switching to flat shading" << std::endl;
				shadingCode = 1;
			}
			else if (event.type == SDL_KEYDOWN &&  event.key.keysym.sym == SDLK_g && shadingCode != 2) {
				std::cout << "switching to gouraud shading" << std::endl;
				shadingCode = 2;
			}
			else if (event.type == SDL_KEYDOWN &&  event.key.keysym.sym == SDLK_p && shadingCode != 3) {
				std::cout << "switching to phong shading" << std::endl;
				shadingCode = 3;
			}
			// Shadow Switch
			if (event.type == SDL_KEYDOWN && event.key.keysym.sym == SDLK_h && shadowCode != 0) {
				std::cout << "turning off shadows" << std::endl;
				shadowCode = 0;
			}
			else if (event.type == SDL_KEYDOWN && event.key.keysym.sym == SDLK_n && shadowCode != 1) {
				std::cout << "switching to hard shadows" << std::endl;
				shadowCode = 1;
			}
			else if (event.type == SDL_KEYDOWN && event.key.keysym.sym == SDLK_m && shadowCode != 2) {
				std::cout << "switching to soft shadows" << std::endl;
				shadowCode = 2;
			}
			// Scene Switch
			if (event.type == SDL_KEYDOWN && event.key.keysym.sym == SDLK_1 && sceneCode != 1) {
				std::cout << "switching to scene 1" << std::endl;
				sceneCode = 1;
			}
			else if (event.type == SDL_KEYDOWN && event.key.keysym.sym == SDLK_2 && sceneCode != 2) {
				std::cout << "switching to scene 2" << std::endl;
				sceneCode = 2;
			}
			else if (event.type == SDL_KEYDOWN && event.key.keysym.sym == SDLK_3 && sceneCode != 3) {
				std::cout << "switching to scene 3" << std::endl;
				sceneCode = 3;
			}
			else if (event.type == SDL_KEYDOWN && event.key.keysym.sym == SDLK_4 && sceneCode != 4) {
				std::cout << "switching to scene 4" << std::endl;
				sceneCode = 4;
			}
			else if (event.type == SDL_KEYDOWN && event.key.keysym.sym == SDLK_5 && sceneCode != 5) {
				std::cout << "switching to scene 5" << std::endl;
				sceneCode = 5;
			}
			if (event.type == SDL_KEYDOWN && event.key.keysym.sym == SDLK_0) {
				std::cout << "starting automatic viewing" << std::endl;
				autoViewing(models, textures, window,focalLength, imageScale, event);
			}

			cameraPosition = handleEvent(event, window, cameraPosition, cameraOrientation);
			cameraOrientation = handleEvent2(event, window, cameraOrientation);
			lightSource = handleEvent3(event, window, lightSource);

			while (event.type == SDL_KEYDOWN && event.key.keysym.sym == SDLK_q ) {
				cameraPosition = orbit(cameraPosition, .005);
				cameraOrientation = lookAt(glm::vec3(0, 0, 0), cameraPosition, cameraOrientation);

				if (renderingCode == 1) wireFrameRenderer(window, cameraPosition, cameraOrientation, focalLength, imageScale, models[sceneCode - 1]);
				else if (renderingCode == 2) rasterisedRenderer(window, cameraPosition, cameraOrientation, focalLength, imageScale, models[sceneCode - 1], textures[0]);
				else if (renderingCode == 3) rayTrace(window, cameraPosition, cameraOrientation, lightSource, focalLength, imageScale, shadingCode, shadowCode, models[sceneCode -1]);
				window.renderFrame();
				window.pollForInputEvents(event);
			}
		}
		if (render) {
			//pointCloudRenderer(window, cameraPosition, focalLength, imageScale, models);
			if (renderingCode == 1) wireFrameRenderer(window, cameraPosition, cameraOrientation, focalLength, imageScale, models[sceneCode - 1]);
			else if (renderingCode == 2) rasterisedRenderer(window, cameraPosition, cameraOrientation, focalLength, imageScale, models[sceneCode - 1], textures[0]);
			else if (renderingCode == 3) rayTrace(window, cameraPosition, cameraOrientation, lightSource, focalLength, imageScale, shadingCode, shadowCode, models[sceneCode -1]);
			window.renderFrame(); // Need to render the frame at the end, or nothing actually gets shown on the screen !
			std::cout << "rendered" << std::endl;
			render = false;
		}
	}
}
