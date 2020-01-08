using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class BSPTree : MonoBehaviour {

    private const float widthSize = 16.0f;
    private const float heightSize = 8.0f;
    private const float maxSize = 2.0f;
    private const float minSize = 0.5f;

    private class Node {
        public Vector2 pos;
        public float width;
        public float height;
        public Node leftChildNode;
        public Node rightChildNode;
        public Node parentNode;
        public bool isConnected;

        public Node(Vector2 pos, float width, float height, Node leftChildNode, Node rightChildNode, Node parentNode) {
            this.pos = pos;
            this.width = width;
            this.height = height;
            this.leftChildNode = leftChildNode;
            this.rightChildNode = rightChildNode;
            this.parentNode = parentNode;
            isConnected = false;
        }

        public void Split() {
            if (width * 0.5f < (maxSize + minSize) / 2 || height * 0.5f < (maxSize + minSize) / 2)
               return;
            //if (width < maxSize && height < maxSize)
             //   return;

            // 判断切割方向
            bool splitVertically;
            if (width > height && width / height > 1.25f) {
                // 纵向切割
                splitVertically = true;
            } else if (height > width && height / width > 1.25f) {
                // 横向切割
                splitVertically = false;
            } else {
                // 随机切割
                splitVertically = Random.Range(0, 100) < 50 ? true : false;
            }
            // 切割
            if (splitVertically) {
                float width1, width2;
                width1 = Random.Range((width / 4), width / 2);
                width2 = Random.Range((width / 4), width / 2);
                float height1, height2;
                height1 = height;
                height2 = height;
                Vector2 pos1 = pos;
                Vector2 pos2 = pos + new Vector2(width / 2, 0);
                leftChildNode = new Node(pos1, width1, height1, null, null, this);
                rightChildNode = new Node(pos2, width2, height2, null, null, this);
            } else {
                float width1, width2;
                width1 = width;
                width2 = width;
                float height1, height2;
                height1 = Random.Range((height / 4), height / 2);
                height2 = Random.Range((height / 4), height / 2);
                Vector2 pos1 = pos;
                Vector2 pos2 = pos + new Vector2(0, height / 2);
                leftChildNode = new Node(pos1, width1, height1, null, null, this);
                rightChildNode = new Node(pos2, width2, height2, null, null, this);
            }
            leftChildNode.Split();
            rightChildNode.Split();
        }

        public void InstantiateRoom(GameObject room, Transform root) {
            if (leftChildNode != null && rightChildNode != null) {
                leftChildNode.InstantiateRoom(room, root);
                rightChildNode.InstantiateRoom(room, root);
            } else if (leftChildNode != null) {
                leftChildNode.InstantiateRoom(room, root);
            } else if (rightChildNode != null) {
                rightChildNode.InstantiateRoom(room, root);
            } else {
                float width = Random.Range(this.width / 2, this.width);
                float height = Random.Range(this.height / 2, this.height);
                Vector2 pos = this.pos + new Vector2(Random.Range(0.0f, this.width - width), Random.Range(0.0f, this.height - height));
                this.pos = pos;
                this.width = width;
                this.height = height;
                room.transform.localScale = new Vector2(width, height);
                Instantiate(room, pos, Quaternion.identity, root);
            }
        }

        public void ConnectRooms(GameObject path, Transform parent) {
            Node anotherNode = GetAnotherRoom();
            if (anotherNode == null) {
                return;
            } else {
                // 连接
                Vector2 pos1 = new Vector2(Random.Range(pos.x, pos.x + width), Random.Range(pos.y, pos.y + height));
                Vector2 pos2 = new Vector2(Random.Range(anotherNode.pos.x, anotherNode.pos.x + anotherNode.width), Random.Range(anotherNode.pos.y, anotherNode.pos.y + anotherNode.height));
                float pathWidth = pos2.x - pos1.x;
                float pathHeight = pos2.y - pos1.y;
                path.transform.localScale = new Vector2(pathWidth, pathHeight);
                Instantiate(path, pos1, Quaternion.identity, parent);
                anotherNode.ConnectRooms(path, parent);
            }
        }

        public Node GetRoom() {
            isConnected = true;
            if (leftChildNode == null && rightChildNode == null) {
                return this;
            } else if (leftChildNode == null) {
                return rightChildNode.GetRoom();
            } else if (rightChildNode == null) {
                return leftChildNode.GetRoom();
            }
            return leftChildNode.GetRoom();
        }

        public Node GetAnotherRoom() {
            if (parentNode == null) {
                if (isConnected && leftChildNode.isConnected && rightChildNode.isConnected) {
                    return null;
                }
            }
            bool atLeftChildNode = this == parentNode?.leftChildNode ? true : false;
            if (atLeftChildNode) {
                if (parentNode?.rightChildNode != null) {
                    if (parentNode.rightChildNode.isConnected) {
                        return parentNode?.GetAnotherRoom();
                    } else {
                        return parentNode?.rightChildNode.GetRoom();
                    }
                } else {
                    return parentNode?.GetAnotherRoom();
                }
            } else {
                if (parentNode?.leftChildNode != null) {
                    if (parentNode.leftChildNode.isConnected) {
                        return parentNode?.GetAnotherRoom();
                    } else {
                        return parentNode?.leftChildNode.GetRoom();
                    }
                } else {
                    return parentNode?.GetAnotherRoom();
                }
            }
        }
    }
    private GameObject room;
    private void Awake() {
        room = Resources.Load<GameObject>("Room");
        Node root = new Node(transform.position, widthSize, heightSize, null, null, null);
        root.Split();
        root.InstantiateRoom(room, transform);
        root.GetRoom().ConnectRooms(room, transform);
    }
}
